#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <conio.h>

#include "fluid_system_host.cuh"		
#include "fluid_system_kern.cuh"
#include "thrust\device_vector.h"	//thrust libs
#include "thrust\sort.h"
#include "fluidMath.cuh"


extern __device__ FluidParams	simData;
extern __device__ ParamCarrier paramCarrier;

//compute drift velocity
__device__ void contributeDriftVel(int i, int muli, cfloat3 ipos, float idens, float ipress, int cell, bufList buf, float* ialpha_pre, float* imassconcen, cfloat3* idriftvelterm, float relax_coef, cfloat3*ialphagrad) {
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist;
	float cmterm;
	float pmterm;
	int j, mulj;

	if (buf.mgridcnt[cell] == 0) return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	cfloat3 force = make_cfloat3(0, 0, 0);
	cfloat3 pgrad[MAX_FLUIDNUM];
	cfloat3 pgradsum;

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];

		if (buf.misbound[j] ==1)
			continue;

		if (buf.MFtype[j] == 1)
			continue;

		mulj = j * MAX_FLUIDNUM;
		dist = (ipos - buf.mpos[j]);		// dist in cm
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale;
		if (dsq < r2 && dsq > 0) {
			//cx = (r2-dsq)*d2;
			dsq = sqrt(dsq*d2);
			c = (simData.psmoothradius - dsq);
			cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j] * buf.mdensity[j];
			//pressure
			pgradsum = make_cfloat3(0, 0, 0);
			for (uint fcount = 0; fcount<simData.mf_catnum; fcount++)
			{
				float jalphaprecount = buf.mf_alpha_pre[mulj+fcount];
				//float ialphaprecount = ialpha_pre[fcount];
				pmterm = cmterm * (-ialpha_pre[fcount]*ipress + jalphaprecount*buf.mpress[j]);
				//pmterm = cmterm * (-ialpha_pre[fcount]*ipress + buf.mf_alpha_pre[mulj+fcount]*buf.mpress[j]);
				pgrad[fcount] = pmterm * dist;
				pgradsum += pgrad[fcount] * imassconcen[fcount];
				//grad alpha
				ialphagrad[fcount] += (jalphaprecount-ialpha_pre[fcount]) * cmterm * dist;
			}

			for (uint fcount = 0; fcount<simData.mf_catnum; fcount++)
			{
				idriftvelterm[fcount] -= relax_coef * (pgrad[fcount]-pgradsum);
			}
		}
	}
}

__global__ void mfComputeDriftVel(bufList buf, int pnum)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= pnum)
		return;

	if (buf.misbound[i] == 1)
		return;
	if (buf.MFtype[i] == 1)
		return;
	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;
	register float relax_coef = simData.relax;					// temporary relax time related coefficient
	register float sigma = simData.mf_diffusion;//0.001f;						//diffusion&tension coefficient
	register float cont, conts, contr;
	cont = simData.cont;
	conts = simData.cont1;
	contr = simData.cont2;
	register cfloat3 accel = buf.mforce[i];				// final accel (g-a) of last step was stored in here cf. advance, 
														//register float massFrack[MAX_FLUIDNUM];
	register uint muloffseti = i * MAX_FLUIDNUM;
	register float invdens = 1.0/buf.mf_restdensity[i];
	register float dsum;
	register float vrx, vry, vrz;
	register float tdiff;
	register cfloat3 ssum;

	register float alpha_pre[MAX_FLUIDNUM], mass_concen[MAX_FLUIDNUM];
	register float ipress = buf.mpress[i];
	register cfloat3 ipos = buf.mpos[i];
	register float idens = buf.mdensity[i];
	register cfloat3 driftVelterm[MAX_FLUIDNUM], alphaGradterm[MAX_FLUIDNUM];
	register cfloat3 sterm[MAX_FLUIDNUM];

	//various viscosity
	relax_coef /= buf.mf_visc[i];

	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{

		alpha_pre[fcount] = buf.mf_alpha_pre[muloffseti+fcount];
		mass_concen[fcount] = alpha_pre[fcount]*simData.mf_dens[fcount]*invdens;
		driftVelterm[fcount] = make_cfloat3(0, 0, 0);
		alphaGradterm[fcount] = make_cfloat3(0, 0, 0);
	}

	for (int c=0; c < simData.gridAdjCnt; c++) {
		contributeDriftVel(i, muloffseti, ipos, idens, ipress, gc + simData.gridAdj[c], buf, alpha_pre, mass_concen, driftVelterm, relax_coef, alphaGradterm);
	}

	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		cfloat3 vel = cont * contr * driftVelterm[fcount];
		buf.mf_vel_phrel[muloffseti+fcount] = vel;

	}

	//first term & second term
	dsum = 0;
	ssum = make_cfloat3(0, 0, 0);
	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		float temp = buf.mf_alpha_pre[muloffseti+fcount];
		dsum += temp * simData.mf_dens[fcount] * simData.mf_dens[fcount] * invdens;

		if (temp>0.0001)
			//sterm[fcount] = buf.mf_alphagrad[muloffseti+fcount]/temp;
			sterm[fcount] = alphaGradterm[fcount]/temp;
		else
			sterm[fcount] = make_cfloat3(0, 0, 0);
		ssum += sterm[fcount] * temp * simData.mf_dens[fcount] * invdens;
	}
	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		tdiff = simData.mf_dens[fcount]-dsum;
		tdiff *= relax_coef;
		vrx = accel.x * tdiff;
		vry = accel.y * tdiff;
		vrz = accel.z * tdiff;
		buf.mf_vel_phrel[muloffseti+fcount] += make_cfloat3(vrx, vry, vrz);

		buf.mf_vel_phrel[muloffseti+fcount] -= cont * conts * sigma * (sterm[fcount]-ssum);
	}


}

//compute alpha advance
__device__ void contributeAlphaChange(int i, int muli, cfloat3 ipos, cfloat3 iveleval, float ipress, float idens, int cell, bufList buf, float* ialpha_pre, float* ialphachange, cfloat3* ivmk)
{
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist, vmr, vkr;
	float cmterm;
	int j, mulj;
	float jalpha_prek;

	if (buf.mgridcnt[cell] == 0) return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];

		if (buf.misbound[j] ==1)
			continue;
		if (buf.MFtype[j] == 1)
			continue;

		mulj = j * MAX_FLUIDNUM;
		dist = (ipos - buf.mpos[j]);		// dist in cm
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale;
		if (dsq < r2 && dsq > 0) {
			dsq = sqrt(dsq*d2);
			c = (simData.psmoothradius - dsq);
			cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j] * buf.mdensity[j];
			vmr = buf.mveleval[j] - iveleval;

			for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
			{
				jalpha_prek = buf.mf_alpha_pre[mulj+fcount];

				ialphachange[fcount] -= 0.5 * cmterm * (jalpha_prek+ialpha_pre[fcount]) * (vmr.x * dist.x + vmr.y * dist.y + vmr.z * dist.z);

				vkr = make_cfloat3((jalpha_prek * buf.mf_vel_phrel[mulj+fcount].x + ialpha_pre[fcount] * ivmk[fcount].x),
					(jalpha_prek * buf.mf_vel_phrel[mulj+fcount].y + ialpha_pre[fcount] * ivmk[fcount].y),
					(jalpha_prek * buf.mf_vel_phrel[mulj+fcount].z + ialpha_pre[fcount] * ivmk[fcount].z));
				ialphachange[fcount] -= cmterm * (vkr.x * dist.x + vkr.y * dist.y + vkr.z * dist.z);
			}
		}
	}
}

__global__ void mfComputeAlphaAdvance(bufList buf, int pnum)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= pnum)
		return;
	if (buf.misbound[i]==1)
		return;
	if (buf.MFtype[i] == 1)
		return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;

	register uint muloffseti = i * MAX_FLUIDNUM;
	register cfloat3 ipos = buf.mpos[i];
	register cfloat3 iveleval = buf.mveleval[i];
	register float ipress = buf.mpress[i];
	register float idens = buf.mdensity[i];
	register float alpha_pre[MAX_FLUIDNUM], alphachange[MAX_FLUIDNUM];
	register cfloat3 ivmk[MAX_FLUIDNUM];

	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		alpha_pre[fcount] = buf.mf_alpha_pre[muloffseti+fcount];
		alphachange[fcount] = 0.0f;
		ivmk[fcount] = buf.mf_vel_phrel[muloffseti+fcount];
	}

	for (int c=0; c < simData.gridAdjCnt; c++) {
		contributeAlphaChange(i, muloffseti, ipos, iveleval, ipress, idens, gc + simData.gridAdj[c], buf, alpha_pre, alphachange, ivmk);
	}

	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		alphachange[fcount] *= simData.mf_dt;

		//alphachange limit
		if (alphachange[fcount]<-0.99)
		{
			alphachange[fcount] = -0.99;// * ((int)(buf.mf_alpha[muloffseti+fcount]>0)-(int)(buf.mf_alpha[muloffseti+fcount]<0));
		}
		buf.mf_alpha[muloffseti+fcount] = alphachange[fcount] + alpha_pre[fcount];

	}



	if (buf.mf_alpha[i*MAX_FLUIDNUM + 2]> simData.tohydro && buf.MFtype[i]==0) { //when sand contains much water
		buf.MFtype[i]=3;	//become fluid particle
		printf("%d become water\n", i);
	}
}

//compute correction
__global__ void mfComputeCorrection(bufList buf, int pnum)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= pnum)
		return;
#ifdef NEW_BOUND
	if (buf.misbound[i]==1)
		return;
#endif

	if (buf.MFtype[i] == 1)
		return;
	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;

	register uint muloffseti = i * MAX_FLUIDNUM;

	float sum;
	float alpha_pre[MAX_FLUIDNUM];
	float alpha_modify;
	int flag;

	sum = 0.0f;
	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		float *temp = buf.mf_alpha + (muloffseti + fcount);
		alpha_pre[fcount] = *temp;
		if (*temp < 0.0000)
			*temp = 0.0f;
		sum += *temp;
	}

	buf.mf_pressure_modify[i] = 0.0f;
	flag = (sum>0.0f);
	sum = flag*sum + (1-flag)*1.0f;

	//int cat = findMaxCat(alpha_pre, simData.mf_catnum, idx, idxlen);
	int maxcat = 0;
	for (uint fcount = 1; fcount<simData.mf_catnum; fcount++)
	{
		if (buf.mf_alpha_pre[i*MAX_FLUIDNUM+fcount]>buf.mf_alpha_pre[i*MAX_FLUIDNUM+maxcat])
			maxcat = fcount;
	}

	sum = 1.0f/sum;
	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		//Correct Alpha
		float* temp = buf.mf_alpha + (muloffseti + fcount);
		*temp = (flag)*(*temp)*sum + (1-flag)*(fcount==maxcat ? 1 : 0);

		//Correct Pressure
		alpha_modify = *temp-alpha_pre[fcount];
		buf.mf_pressure_modify[i] += -(6*pow(1/(buf.mdensity[i]*buf.mf_restdensity[i]), 7)+1.0) * simData.pintstiff * simData.mf_dens[fcount] * alpha_modify;
	}

	//Update Color
	//buf.mclr[i] = COLORA( buf.mf_alpha[i*MAX_FLUIDNUM],  buf.mf_alpha[i*MAX_FLUIDNUM+1],  buf.mf_alpha[i*MAX_FLUIDNUM+2], 1);

	////Restdensity Update
	float newdens = 0.0;
	float newvisc = 0.0;
	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		newdens += buf.mf_alpha[i*MAX_FLUIDNUM+fcount] * simData.mf_dens[fcount];
		newvisc += buf.mf_alpha[i*MAX_FLUIDNUM+fcount] * simData.mf_visc[fcount];
	}

	if (buf.MFtype[i]!=2 && buf.MFtype[i]!=1)
		buf.mf_restdensity[i] = newdens;
	buf.mf_visc[i] = newvisc;
}

__device__ float gamma(float q)
{
	if (q<2.0/3.0 && q>0)
		return 2.0/3.0;
	if (q>=2.0/3.0 && q<1)
		return 2*q-3.0/2.0*q*q;
	if (q>=1 && q<2)
		return (2-q)*(2-q)/2.0;
	return 0;
}

//advance particles
__device__ void mfChRebalance(int i, int muli, bufList buf, int firstReactor, int secondReactor, int product)
{
	float chGamma = 0.01;
	register float alpha1 = buf.mf_alpha[muli+firstReactor];
	register float alpha2 = buf.mf_alpha[muli+secondReactor];
	//register float alphap;
	register float massTrans1, massTrans2;
	//register float V0 = buf.mf_restmass[i] * buf.mdensity[i];
	register float Vp;
	register float rhop1 = simData.mf_dens[firstReactor];
	register float rhop2 = simData.mf_dens[secondReactor];
	register float rhopp = simData.mf_dens[product];
	register float deltaAlphaP;

	//chGamma *= (alpha1*alpha2);
	chGamma *= (alpha1+alpha2);
	if (chGamma == 0)return;
	if (chGamma > alpha1)chGamma = alpha1;
	if (chGamma > alpha2)chGamma = alpha2;

	massTrans1 = chGamma * rhop1;
	massTrans2 = chGamma * rhop2;

	deltaAlphaP = (massTrans1 + massTrans2) / rhopp;

	Vp = 1 + deltaAlphaP - 2 * chGamma;
	Vp = 1/Vp;
	buf.mf_alpha[muli+firstReactor] -= chGamma;
	buf.mf_alpha[muli+secondReactor] -= chGamma;
	buf.mf_alpha[muli+product] += deltaAlphaP;

	for (uint fcount = 0; fcount<simData.mf_catnum; fcount++)
	{
		buf.mf_alpha[muli+fcount] *= Vp;
	}

	buf.mf_restdensity[i] *= Vp;
}


__device__ cfloat3 contributeForce(int i, int muli, cfloat3 ipos, cfloat3 iveleval, float ipress, float idens, int cell, bufList buf, cfloat3* ivelxcor, float ivisc)
{
	//Force here represents the acceleration
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist, vmr;
	float cmterm, cmterm1;

	float pmterm, vmterm;

	int j;
	float aveDenij, cx, xterm;

	if (buf.mgridcnt[cell] == 0) return make_cfloat3(0, 0, 0);

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	cfloat3 force = make_cfloat3(0, 0, 0);

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];
		//if (buf.misbound[j]==1)
		//    continue;

		dist = (ipos - buf.mpos[j]);		// dist in cm
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale;

		if (dsq < r2 && dsq > 0) {
			cx = (r2-dsq)*d2;
			dsq = sqrt(dsq*d2);
			c = (simData.psmoothradius - dsq);

			cmterm1 = simData.spikykern * c * c / dsq;
			cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j] * buf.mdensity[j];


			if (buf.misbound[j] != 1) //fluid & solid, fluid & fluid, solid & solid
			{
				//pressure
				pmterm = -0.5f * cmterm * (ipress  + buf.mpress[j])*idens;
				force += pmterm * dist;

				//viscosity
				vmr = iveleval - buf.mveleval[j];
				vmterm = cmterm * (ivisc+buf.mf_visc[j]) * idens;
				force += vmterm * vmr;

			}
			else { //boundary & non-boundary
				   //pressure
				   //pmterm = - cmterm1 * buf.mf_restdensity[i] * buf.mf_restmass[j] / buf.mf_restdensity[j] * (ipress) *buf.mdensity[i]*buf.mdensity[i];

				pmterm = -0.5f * cmterm * (ipress  + buf.mpress[j])*idens;

				force += pmterm * dist;

				//viscosity prevent penetration
				vmr = iveleval - buf.mveleval[j];
				float pi_ij = vmr.x*dist.x + vmr.y*dist.y + vmr.z*dist.z;
				if (pi_ij < 0) {
					pi_ij = pi_ij / (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z + r2 * 0.01);
					pi_ij = pi_ij * 2 * simData.psmoothradius * (ivisc + buf.mf_visc[j]) * idens /2;
					pi_ij = - cmterm1 * buf.mf_restdensity[i] * buf.mf_restmass[j]/buf.mf_restdensity[j] * pi_ij;
					force += pi_ij * dist * simData.visc_factor;
				}
			}

			//XSPH correction
			aveDenij = 2/(1/buf.mdensity[j]+1/idens);
			xterm = cx*cx*cx*buf.mf_restmass[j]*aveDenij*simData.poly6kern*0.5; //0.5=epsilon
			ivelxcor->x += -vmr.x * xterm;
			ivelxcor->y += -vmr.y * xterm;
			ivelxcor->z += -vmr.z * xterm;
		}
	}
	return force;
}



// **********   Project-u  Functions *********
__device__ cfloat3 contributeForce_projectu(int i, int muli, cfloat3 ipos, cfloat3 iveleval, float ipress, float idens, int cell, bufList buf, float* ialpha_pre, float ipressure_modify, cfloat3* ivmk, cfloat3* ivelxcor, float ivisc)
{
	//Force here represents the acceleration
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist, vmr;
	float cmterm, cmterm1;

	float pmterm, vmterm;

	int j, mulj;
	float aveDenij, cx, xterm;

	if (buf.mgridcnt[cell] == 0) return make_cfloat3(0, 0, 0);

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	cfloat3 force = make_cfloat3(0, 0, 0);

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];

		mulj = j * MAX_FLUIDNUM;
		dist = (ipos - buf.mpos[j]);		// dist in cm
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale;

		if (dsq < r2 && dsq > 0) {
			cx = (r2-dsq)*d2;
			dsq = sqrt(dsq*d2);
			c = (simData.psmoothradius - dsq);

			cmterm1 = simData.spikykern * c * c / dsq;
			cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j] * buf.mdensity[j];
			//pressure

			if (buf.misbound[j] != 1) //fluid & solid, fluid & fluid, solid & solid
			{
				pmterm = -0.5f * cmterm * (ipress  + buf.mpress[j])*idens;
				force += pmterm * dist;

				//viscosity
				vmr = iveleval - buf.mveleval[j];
				vmterm = cmterm * (ivisc+buf.mf_visc[j]) * idens;
				force += vmterm * vmr;

			}
			else { //boundary & non-boundary

				   //pressure
				pmterm = - cmterm1 * buf.mf_restdensity[i] * buf.mf_restmass[j] / buf.mf_restdensity[j] * (ipress)*buf.mdensity[i]*buf.mdensity[i];
				force += pmterm * dist * simData.omega;

				//viscosity
				vmr = iveleval - buf.mveleval[j];
				float pi_ij = vmr.x*dist.x + vmr.y*dist.y + vmr.z*dist.z;
				if (pi_ij < 0) {
					pi_ij = pi_ij / (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z + r2 * 0.01);
					pi_ij = pi_ij * 2 * simData.psmoothradius * (ivisc + buf.mf_visc[j]) * idens /2;
					pi_ij = - cmterm1 * buf.mf_restdensity[i] * buf.mf_restmass[j]/buf.mf_restdensity[j] * pi_ij;
					force += pi_ij * dist * simData.visc_factor;
				}
			}

			if (buf.misbound[j] != 1) //drift velocity tensor
									  //T_dm
				for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
				{
					cfloat3 dtermj = cmterm * (buf.mf_vel_phrel[mulj+fcount].x * dist.x + buf.mf_vel_phrel[mulj+fcount].y * dist.y + buf.mf_vel_phrel[mulj+fcount].z * dist.z) * buf.mf_alpha_pre[mulj+fcount] * buf.mf_vel_phrel[mulj+fcount];
					cfloat3 dtermi = cmterm * (ivmk[fcount].x * dist.x + ivmk[fcount].y * dist.y + ivmk[fcount].z * dist.z) * ialpha_pre[fcount] * ivmk[fcount];
					force += (dtermj + dtermi) * simData.mf_dens[fcount] * idens;
				}

#ifndef _nXSPH
			//XSPH correction
			aveDenij = 2/(1/buf.mdensity[j]+1/idens);
			xterm = cx*cx*cx*buf.mf_restmass[j]*aveDenij*simData.poly6kern*0.5; //0.5=epsilon
			ivelxcor->x += -vmr.x * xterm;
			ivelxcor->y += -vmr.y * xterm;
			ivelxcor->z += -vmr.z * xterm;
		}
#endif

	}
	return force;
}

__global__ void ComputeForce_projectu(bufList buf, int pnum)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= pnum)
		return;
	if (buf.misbound[i]==1)
		return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;

	register uint muloffseti = i * MAX_FLUIDNUM;
	register cfloat3 ipos = buf.mpos[i];
	register cfloat3 iveleval = buf.mveleval[i];
	register float ipress = buf.mpress[i];
	register float idens = buf.mdensity[i];
	register float alpha_pre[MAX_FLUIDNUM];
	register cfloat3 ivmk[MAX_FLUIDNUM];
	register float pressure_modify = buf.mf_pressure_modify[i];
	register cfloat3 *ivelxcor = buf.mf_velxcor+i;
	register float ivisc = buf.mf_visc[i];

	register cfloat3 force = make_cfloat3(0, 0, 0);
	*ivelxcor = make_cfloat3(0, 0, 0);

	for (uint fcount = 0; fcount < simData.mf_catnum; fcount++)
	{
		alpha_pre[fcount] = buf.mf_alpha_pre[muloffseti+fcount];
		ivmk[fcount] = buf.mf_vel_phrel[muloffseti+fcount];
	}

	for (int c=0; c < simData.gridAdjCnt; c++) {
		force += contributeForce_projectu(i, muloffseti, ipos, iveleval, ipress, idens, gc + simData.gridAdj[c], buf, alpha_pre, pressure_modify, ivmk, ivelxcor, ivisc);
	}
	buf.mforce[i] = force;
}

__device__ void contributeVelocityGradient(float* result, int i, cfloat3 ipos, cfloat3 iveleval, int cell, bufList buf)
{
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist, jveleval;
	float cmterm;
	//	float massj,massi;

	//	float q;
	int j;
	//	float aveDenij,cx,xterm;
	float alphaij;
	if (buf.mgridcnt[cell] == 0) return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	//massi = buf.mf_restmass[i];
	int mftype = buf.MFtype[i];

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];

		//massj = buf.mf_restmass[j];
		//jveleval = buf.mveleval[j]*buf.mdensity[j]*buf.mdensity[j] + iveleval*buf.mdensity[i]*buf.mdensity[i];
		jveleval = buf.mveleval[j]+buf.mf_vel_phrel[j*MAX_FLUIDNUM] - iveleval-buf.mf_vel_phrel[i*MAX_FLUIDNUM];

		//alphaij = sqrt(buf.mf_alpha[i*MAX_FLUIDNUM + mftype] * buf.mf_alpha[j*MAX_FLUIDNUM + mftype]);

		if (buf.mf_alpha[i*MAX_FLUIDNUM + mftype]<0.000001 || buf.mf_alpha[j*MAX_FLUIDNUM + mftype]<0.000001)
			alphaij = 0;
		else
			alphaij = 2.0 / (1.0/buf.mf_alpha[i*MAX_FLUIDNUM + mftype] + 1.0/buf.mf_alpha[j*MAX_FLUIDNUM+mftype]);

		dist = (ipos - buf.mpos[j]);		// dist in cm
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale*alphaij;

		if (dsq < r2 && dsq > 0) {
			dsq = sqrt(dsq * d2);
			c = (simData.psmoothradius - dsq);
			cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j] * buf.mdensity[j];
			//cmterm = simData.spikykern * c * c / dsq;
			jveleval = jveleval * cmterm;
			result[0] += jveleval.x * dist.x;	result[1] += jveleval.x * dist.y;	result[2] += jveleval.x * dist.z;
			result[3] += jveleval.y * dist.x;	result[4] += jveleval.y * dist.y;	result[5] += jveleval.y * dist.z;
			result[6] += jveleval.z * dist.x;	result[7] += jveleval.z * dist.y;	result[8] += jveleval.z * dist.z;
		}
	}
}

__device__ void print9(char* string, int id, float* buf) {

	printf("%s %d \n%f %f %f\n%f %f %f\n%f %f %f\n", string, id, buf[0], buf[1], buf[2],
		buf[3], buf[4], buf[5], buf[6], buf[7], buf[8]);
	return;
}

__global__ void ComputeSPHtensor(bufList buf, int pnum) { //update sigma
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= pnum)
		return;

	if (buf.misbound[i]==1) { //boundary
		float* sigma = buf.MFtemptensor+i*9;
		float* sigma_s = buf.MFtensor + i*9;
		for (int k=0; k<9; k++) {
			sigma[k] = 0;
			sigma_s[k] = 0;
		}
		sigma[0] -= buf.mpress[i];// +buf.mf_pressure_modify[i];
		sigma[4] -= buf.mpress[i];// +buf.mf_pressure_modify[i];
		sigma[8] -= buf.mpress[i];// +buf.mf_pressure_modify[i];
		float* rtensor = buf.MFRtensor+i*9;
		for (int k=0; k<9; k++) {
			rtensor[k] = sigma[k];
			if (rtensor[k]>0)
				rtensor[k]=0;
		}

		return;
	}
	// Get search cell
	int nadj = (1 * simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;

	if (buf.MFtype[i] == 0 || buf.MFtype[i]==3) // let 3 be water
	{
		register cfloat3 ipos = buf.mpos[i], iveleval = buf.mveleval[i];
		float epsilon[9], omega[9], vgrad[9], vgradtrans[9];
		float dsigma[9];
		register float lambda, trepsilon, sqrtj2;

		float coG = simData.coG;//constants
		float phi = simData.phi;

		//calculate velocity gradient
		float* mfvgrad = buf.MFvelgrad + i * 9;
		for (int k = 0; k<9; k++)
			mfvgrad[k] = 0;
		for (int c = 0; c < simData.gridAdjCnt; c++)
			contributeVelocityGradient(mfvgrad, i, ipos, iveleval, gc + simData.gridAdj[c], buf);
		for (int k = 0; k<9; k++)
			vgrad[k] = mfvgrad[k];

		//calculate epsilon and omega
		transpose3(vgrad, vgradtrans);
		for (int k = 0; k<9; k++) {
			epsilon[k] = (vgrad[k] + vgradtrans[k])*0.5;
			omega[k] = (vgrad[k] - vgradtrans[k])*0.5;
		}
		trepsilon = epsilon[0] + epsilon[4] + epsilon[8];
		buf.MFpepsilon[i] = trepsilon / 3.f;
		float* sigma_s = buf.MFtensor + i * 9;
		float* sigma = buf.MFtemptensor + i * 9;

		//get sigma again
		for (int k = 0; k<9; k++) {
			sigma[k] = sigma_s[k];
		}
		sigma[0] -= buf.mpress[i];
		sigma[4] -= buf.mpress[i];
		sigma[8] -= buf.mpress[i];

		//get lambda
		float u[9], v[9], vt[9], s[9], comp[3];

		//calculate e_ab
		trepsilon = (epsilon[0] + epsilon[4] + epsilon[8]);
		for (int k = 0; k<9; k++)
			vgradtrans[k] = epsilon[k];
		vgradtrans[0] -= trepsilon / 3;
		vgradtrans[4] -= trepsilon / 3;
		vgradtrans[8] -= trepsilon / 3;

		float tgphi, yieldf;
		float coLambdaK = simData.coLambdaK, alpha_phi, kc, cohesion = simData.cohesion;
		tgphi = tan(phi / 180 * 3.14159265359);
		alpha_phi = phi;// tgphi / (sqrtf(9 + 12 * tgphi*tgphi));
		kc = 3 * cohesion / sqrtf(9 + 12 * tgphi*tgphi);
		//calculate lambda
		sqrtj2 = 0;
		lambda = 0;
		for (int k = 0; k<9; k++) {
			sqrtj2 += sigma_s[k] * sigma_s[k];
		}
		sqrtj2 = sqrtf(sqrtj2*0.5);
		buf.mclr[i] = COLORA(1, 1, 1, 1);

		if (alpha_phi*buf.mpress[i] * 3 + kc<sqrtj2)
		{
			buf.mclr[i] = COLORA(1, 1, 1, 1);
			float rn = 0;
			if (sqrtj2>1e-4)
				rn = (alpha_phi*buf.mpress[i] * 3 + kc) / sqrtj2;
			for (int k = 0; k<9; k++)
				sigma_s[k] *= rn;
			for (int k = 0; k<9; k++) {
				//sab_epab += sigma_s[k]*epsilon[k];
				sqrtj2 += sigma_s[k] * sigma_s[k];
				//lambda += sigma_s[k]*vgradtrans[k];
			}
			sqrtj2 = sqrtf(sqrtj2*0.5);
		}
		if (sqrtj2 >1e-4)
		{
			for (int k = 0; k < 9; k++)
			{
				lambda += (coG / sqrtj2)*sigma_s[k] * vgradtrans[k];
			}
			lambda += 9 * alpha_phi*(buf.mpress[i] - buf.last_mpress[i]) / simData.mf_dt;
			lambda /= 9 * alpha_phi*alpha_phi*coLambdaK + coG;
		}
		else lambda = 0;
		lambda *= simData.Yradius;
		buf.last_mpress[i] = buf.mpress[i];

		//update s
		float norm_depsilon = 0;
		for (int k = 0; k < 9; k++)
			norm_depsilon += epsilon[k] * epsilon[k];
		norm_depsilon = sqrtf(norm_depsilon);
		//if (i % 100000000 == 0){
		//	for (int k = 0; k < 9; k++)
		//		printf("%f ", norm_depsilon*(sigma_s[k] + sigma_s[k])*simData.coD);
		//	printf("\n");
		//}
		for (int ii = 0; ii<3; ii++)
			for (int ij = 0; ij<3; ij++) {
				dsigma[ii * 3 + ij] = 0;
				for (int ik = 0; ik<3; ik++) {
					dsigma[ii * 3 + ij] += sigma[ii * 3 + ik] * omega[ij * 3 + ik];
					dsigma[ii * 3 + ij] += sigma[ik * 3 + ij] * omega[ii * 3 + ik];
				}
				dsigma[ii * 3 + ij] += 2 * coG*buf.mpress[i] * vgradtrans[ii * 3 + ij] + norm_depsilon*(sigma_s[ii * 3 + ij] + sigma_s[ii * 3 + ij])*simData.coD;
				if (fabs(sigma[0] + sigma[4] + sigma[8])>1e-4)
				{
					float trdouble = 0;
					for (int k = 0; k < 9; k++)
						trdouble += sigma[k] * epsilon[k];
					dsigma[ii * 3 + ij] += simData.coD0*trdouble / (sigma[0] + sigma[4] + sigma[8]) * sigma[ii * 3 + ij];
				}

				if (fabs(lambda) > 1e-5)
				{
					dsigma[ii * 3 + ij] -= lambda * coG / sqrtj2*sigma_s[ii * 3 + ij];
					if (ii == ij)
						dsigma[ii * 3 + ij] -= lambda*coLambdaK*alpha_phi * 3;
				}

			}
		//printf("%f\n", lambda);
		for (int k = 0; k<9; k++) {
			sigma_s[k] = sigma_s[k] + dsigma[k] * simData.mf_dt;
			sigma[k] = sigma_s[k];

		}
		float tmp = sigma_s[0] + sigma_s[4] + sigma_s[8];
		tmp /= 3;
		sigma_s[0] -= tmp;
		sigma_s[4] -= tmp;
		sigma_s[8] -= tmp;
		//calculate sigma
		sigma[0] -= buf.mpress[i];
		sigma[4] -= buf.mpress[i];
		sigma[8] -= buf.mpress[i];

		//Calculate artificial tensor

		for (int k = 0; k<9; k++) {
			u[k] = sigma[k];
			s[k] = 0;
		}
		svdecomp3(comp, u, v, 0.000001);

		for (int ii = 0; ii<3; ii++) {
			if (u[ii * 3] * v[ii * 3]<0) {
				for (int ij = 0; ij<3; ij++) {
					v[ii * 3 + ij] *= -1;
				}
				comp[ii] *= -1;
			}
		}

		for (int k = 0; k<3; k++) {
			if (comp[k]<0)
				comp[k] = 0;
		}
		s[0] = comp[0]; s[4] = comp[1]; s[8] = comp[2];
		transpose3(v, vt);
		multiply_matrix3(u, s, u);
		multiply_matrix3(u, vt, u);

		float* rtensor = buf.MFRtensor + i * 9;
		for (int k = 0; k<9; k++)
			rtensor[k] = u[k];
		//if (i % 1000 == 0)
		//	printf("shit1");
	}
	else if (buf.MFtype[i] == 1)
	{
		register cfloat3 ipos = buf.mpos[i], iveleval = buf.mveleval[i];
		float epsilon[9], omega[9], vgrad[9], vgradtrans[9];
		float dsigma[9];
		register float lamda, trepsilon, j2, tempf;
		//register float sab_epab, trsigma, alpha_phi,kc;
		float coG = simData.solid_coG;//constants
									  //float phi = simData.phi;

									  //calculate velocity gradient
		float* mfvgrad = buf.MFvelgrad + i * 9;
		for (int k = 0; k<9; k++)
			mfvgrad[k] = 0;
		for (int c = 0; c < simData.gridAdjCnt; c++)
			contributeVelocityGradient(mfvgrad, i, ipos, iveleval, gc + simData.gridAdj[c], buf);
		for (int k = 0; k<9; k++)
			vgrad[k] = mfvgrad[k];

		//calculate epsilon and omega
		transpose3(vgrad, vgradtrans);
		for (int k = 0; k<9; k++) {
			epsilon[k] = (vgrad[k] + vgradtrans[k])*0.5;
			omega[k] = (vgrad[k] - vgradtrans[k])*0.5;
		}
		trepsilon = epsilon[0] + epsilon[4] + epsilon[8];
		buf.MFpepsilon[i] = trepsilon / 3.f;
		float* sigma_s = buf.MFtensor + i * 9;
		float* sigma = buf.MFtemptensor + i * 9;

		//get sigma again
		for (int k = 0; k<9; k++) {
			sigma[k] = sigma_s[k];
		}
		sigma[0] -= buf.mpress[i];
		sigma[4] -= buf.mpress[i];
		sigma[8] -= buf.mpress[i];

		//get lamda
		float u[9], v[9], vt[9], s[9], comp[3];

		//calculate e_ab
		trepsilon = (epsilon[0] + epsilon[4] + epsilon[8]);
		for (int k = 0; k<9; k++)
			vgradtrans[k] = epsilon[k];
		vgradtrans[0] -= trepsilon / 3;
		vgradtrans[4] -= trepsilon / 3;
		vgradtrans[8] -= trepsilon / 3;

		//calculate lamda
		j2 = 0;
		lamda = 0;
		for (int k = 0; k<9; k++) {
			//sab_epab += sigma_s[k]*epsilon[k];
			j2 += 0.5*sigma_s[k] * sigma_s[k];
			lamda += sigma_s[k] * vgradtrans[k];
		}

		if (j2 / 2 / simData.solid_coG - simData.solid_Yradius >0) {
			tempf = simData.solid_Yradius * 2 * simData.solid_coG;
			tempf = sqrt(tempf / j2);
			for (int k = 0; k<9; k++) {
				sigma_s[k] = sigma_s[k] * tempf;
			}
			j2 = simData.solid_Yradius * 2 * simData.solid_coG;
			lamda *= tempf;
		}

		if (abs(lamda)>0.00001)
			lamda /= (simData.solid_Yradius * j2);

		//update s
		for (int ii = 0; ii<3; ii++)
			for (int ij = 0; ij<3; ij++) {
				dsigma[ii * 3 + ij] = 0;
				for (int ik = 0; ik<3; ik++) {
					dsigma[ii * 3 + ij] += sigma[ii * 3 + ik] * omega[ij * 3 + ik];
					dsigma[ii * 3 + ij] += sigma[ik * 3 + ij] * omega[ii * 3 + ik];
				}
				dsigma[ii * 3 + ij] += 2 * coG * vgradtrans[ii * 3 + ij];

				dsigma[ii * 3 + ij] -= lamda * sigma_s[ii * 3 + ij];

			}

		for (int k = 0; k<9; k++) {
			sigma_s[k] = sigma_s[k] + dsigma[k] * simData.mf_dt;
			sigma[k] = sigma_s[k];
			if (fabs(sigma[k])<1e-6)
				sigma[k] = 0;
		}

		//calculate sigma
		sigma[0] -= buf.mpress[i];
		sigma[4] -= buf.mpress[i];
		sigma[8] -= buf.mpress[i];

		//Calculate artificial tensor
		float total = 0;
		for (int k = 0; k<9; k++) {
			u[k] = sigma[k];
			total += fabs(sigma[k]);
			s[k] = 0;
		}
		if (total < 1e-5)
		{
			float* rtensor = buf.MFRtensor + i * 9;
			for (int k = 0; k<9; k++)
				rtensor[k] = u[k];
			return;
		}
		svdecomp3(comp, u, v, 0.000001);

		for (int ii = 0; ii<3; ii++) {
			if ((u[ii] + u[ii + 3] + u[ii + 6]) * (v[ii] + v[ii + 3] + v[ii + 6])<0) {
				for (int ij = 0; ij<3; ij++) {
					v[ij * 3 + ii] *= -1;
				}
				comp[ii] *= -1;
			}
		}

		for (int k = 0; k<3; k++) {
			if (comp[k]<0)
				comp[k] = 0;
		}
		s[0] = comp[0]; s[4] = comp[1]; s[8] = comp[2];
		transpose3(v, vt);
		multiply_matrix3(u, s, u);
		multiply_matrix3(u, vt, u);

		float* rtensor = buf.MFRtensor + i * 9;
		for (int k = 0; k<9; k++)
			rtensor[k] = u[k];
	}
}

__device__ cfloat3 contributeSigmaForce(int i, cfloat3 ipos, float* itensor, int cell, bufList buf)
{
	if (buf.mgridcnt[cell] == 0)
		return make_cfloat3(0, 0, 0);
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist, origin_dist;
	float* jtensor;
	float cmterm;
	//	float massj,massi;
	int j;
	cfloat3 vij, xij;
	float xvprod, densityij, phiij;
	float PIij;
	cfloat3 visc;
	float fij, rij;
	//float Wij,Wdh;

	cfloat3 result = make_cfloat3(0, 0, 0);

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	//massi = buf.mf_restmass[i];
	float tensorij[9];
	for (int k=0; k<9; k++)
		tensorij[k] = 0;
	//float Rij[9];
	float alphaij;
	int mftype = buf.MFtype[i];
	int mftypej;

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];
		mftypej = buf.MFtype[j];

		dist = (ipos - buf.mpos[j]);		// dist in cm
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale;
		origin_dist = -dist;
		if (!(dsq < r2 && dsq > 0))
			continue;

		vij =  buf.mveleval[i] - buf.mveleval[j];
		xij =  (ipos - buf.mpos[j])*simData.psimscale;
		rij = sqrt(xij.x*xij.x + xij.y*xij.y + xij.z*xij.z);
		fij = pow((simData.psmoothradius- rij)/(simData.psmoothradius- simData.initspacing),
			3.0f * simData.coN);
		float rfactor = buf.mf_restdensity[i] - buf.mf_restdensity[j];
		if (rfactor<0) rfactor*= -1;
		rfactor = rfactor * 0.08 / (buf.mf_restdensity[i] + buf.mf_restdensity[j]) * buf.mdensity[i] * buf.mdensity[j];
		float rfactor_inter = 1.0*buf.mdensity[i] * buf.mdensity[j];

		jtensor = buf.MFtensor + j*9;

		if (mftype==0 || mftype==1) {
			if (buf.mf_alpha[i*MAX_FLUIDNUM + mftype]<0.000001 || buf.mf_alpha[j*MAX_FLUIDNUM + mftypej]<0.000001)
				alphaij = 0;
			else
				alphaij = 2.0 / (1.0/buf.mf_alpha[i*MAX_FLUIDNUM + mftype] + 1.0/buf.mf_alpha[j*MAX_FLUIDNUM+mftypej]);
			//alphaij = sqrt(buf.mf_alpha[i*MAX_FLUIDNUM + mftype] * buf.mf_alpha[j*MAX_FLUIDNUM + mftype]);
		}
		else
			alphaij = 0;

		//sigma force related to pressure
		//R term involved to prevent particle clumping (not working)
		for (int k=0; k<9; k++) {
			if (mftype==0 || mftype==1)
				tensorij[k] = (jtensor[k]*buf.mdensity[j]*buf.mdensity[j] + itensor[k]*buf.mdensity[i]*buf.mdensity[i])*alphaij;
			else
				tensorij[k] = 0;

			//Handle density ratio
			//if(buf.MFtype[i] != buf.MFtype[j])
			//	tensorij[k] += fij*(buf.MFRtensor[i*9+k]+buf.MFRtensor[j*9+k])*simData.pan_r * rfactor;//R term

			if (buf.MFtype[i] == 1 && buf.MFtype[j] == 1)//Artificial tensor
				tensorij[k] -= fij*(buf.MFRtensor[i * 9 + k] * buf.mdensity[i] * buf.mdensity[i] + buf.MFRtensor[j * 9 + k] * buf.mdensity[j] * buf.mdensity[j])*simData.pan_r;
		}
		dsq = sqrt(dsq * d2);
		c = (simData.psmoothradius - dsq);
		cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j];
		dist *= cmterm;

		result.x += tensorij[0]*dist.x + tensorij[1]*dist.y + tensorij[2]*dist.z;
		result.y += tensorij[3]*dist.x + tensorij[4]*dist.y + tensorij[5]*dist.z;
		result.z += tensorij[6]*dist.x + tensorij[7]*dist.y + tensorij[8]*dist.z;

		//viscosity force
		xvprod = xij.x * vij.x + xij.y * vij.y + xij.z * vij.z;
		if (xvprod < 0) {
			//if(buf.MFtype[i]==1){
			phiij = simData.psmoothradius * xvprod/
				((xij.x*xij.x + xij.y*xij.y + xij.z*xij.z) + 0.01* simData.r2);
			densityij = (1/buf.mdensity[i]+1/buf.mdensity[j])*0.5;

			//if( buf.MFtype[i]==1 && buf.MFtype[j]==1) //inner visc of sand
			//	PIij = (-simData.coA * phiij + simData.coB * phiij*phiij)/densityij;

			//if ( buf.MFtype[i]!=0 && buf.misbound[j] == 1 || (buf.MFtype[j] != 1 && buf.MFtype[i] == 1) || (buf.MFtype[j] == 1 && buf.MFtype[i] != 1))	//force from boundary to prevent penetrating
			if (buf.misbound[j] == 1 || (buf.MFtype[j] != 1 && buf.MFtype[i] == 1) || (buf.MFtype[j] == 1 && buf.MFtype[i] != 1))	//force from boundary to prevent penetrating	
				PIij = (-simData.fsa * phiij + simData.fsb * phiij*phiij) / densityij;
			else if (buf.MFtype[j] == buf.MFtype[i])
				PIij = (-simData.coA * phiij + simData.coB * phiij*phiij)/densityij * alphaij;
			else if (buf.MFtype[i] == 1 && buf.MFtype[j] == 0 || buf.MFtype[i] == 0 && buf.MFtype[j] == 1)
				PIij = (-simData.fsa * phiij + simData.fsb * phiij*phiij) / densityij;
			else
				PIij = 0;
			if (buf.misbound[j] != 0 && buf.mpos[j].z < 3)
				result = result - buf.mf_alpha[i*MAX_FLUIDNUM] * simData.boundaryVisc*buf.mveleval[i] * dot(origin_dist, dist) / dot(origin_dist, origin_dist)*buf.mdensity[i] * buf.mf_restdensity[i] * buf.mf_restmass[j] / buf.mf_restdensity[j];

			//if( buf.MFtype[j]==2 && buf.MFtype[i]==1) //visc force from boundary to sand
			//	result = result -simData.boundaryVisc*buf.mveleval[i] * dot(origin_dist, dist) / dot(origin_dist, origin_dist)*buf.mdensity[i]* buf.mf_restdensity[i] * buf.mf_restmass[j] / buf.mf_restdensity[j];

			/*else if (buf.MFtype[i] != buf.MFtype[j])
			{
			PIij = (-simData.fsa * phiij + simData.fsb * phiij*phiij) / densityij;
			if (buf.mpos[j].y<5)
			result = result -simData.boundaryVisc*buf.mveleval[i] * dot(origin_dist, dist) / dot(origin_dist, origin_dist)*buf.mdensity[i]* buf.mf_restdensity[i] * buf.mf_restmass[j] / buf.mf_restdensity[j];
			}
			else
			PIij = 0;*/
			result.x -= PIij * dist.x;
			result.y -= PIij * dist.y;
			result.z -= PIij * dist.z;
		}

	}

	return result;
}

__global__ void AddSPHtensorForce(bufList buf, int pnum) {
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= pnum) return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;

	register cfloat3 ipos = buf.mpos[i];
	float *itensor = buf.MFtensor + i*9;
	cfloat3 tensorForce =  make_cfloat3(0, 0, 0);

	for (int c=0; c < simData.gridAdjCnt; c++) {
		tensorForce += contributeSigmaForce(i, ipos, itensor, gc+simData.gridAdj[c], buf);
	}
	buf.mforce[i] = buf.mforce[i] + tensorForce;

}

//**********************  end project-u    ************************

