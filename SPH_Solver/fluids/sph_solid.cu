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


//__device__ float mfContributePressureInit(int i, cfloat3 p, int cell, bufList buf)
//{
//	cfloat3 dist;
//	float dsq, c, sum;
//	float massj;
//	register float d2 = simData.psimscale * simData.psimscale;
//	register float r2 = simData.r2/d2;
//
//	sum = 0.0;
//	int j;
//
//	if (buf.mgridcnt[cell] == 0)
//		return 0.0;
//
//	int cfirst = buf.mgridoff[cell];
//	int clast = cfirst + buf.mgridcnt[cell];
//
//	for (int cndx = cfirst; cndx < clast; cndx++) {
//		j = buf.mgrid[cndx];
//
//		if (buf.MFtype[i]!=buf.MFtype[j])
//			continue;
//
//		dist = p - buf.mpos[j];
//		massj = buf.mf_restmass[j];
//		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
//
//		if (dsq < r2 && dsq > 0.0) {
//			c = (r2 - dsq)*d2;
//			sum += c * c * c * massj;
//		}
//	}
//
//	return sum;
//}
//
//__global__ void initDensity(bufList buf, int pnum) {
//
//	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
//	if (i >= pnum)
//		return;
//
//	if (buf.MFtype[i]!= 1) //solid only
//		return;
//
//	// Get search cell
//	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
//	uint gc = buf.mgcell[i];
//	if (gc == GRID_UNDEF) return;						// particle out-of-range
//	gc -= nadj;
//
//	cfloat3 pos = buf.mpos[i];
//	float sum = 0.0;
//	for (int c=0; c < simData.gridAdjCnt; c++) {
//		sum += mfContributePressureInit(i, pos, gc + simData.gridAdj[c], buf);
//	}
//	sum += simData.r2 * simData.r2 * simData.r2 * buf.mf_restmass[i];
//	sum = sum * simData.poly6kern;
//	//now sum is density
//	buf.mf_restdensity[i] = sum;
//}

//__device__ void contributeVGrad(float* result, int i, cfloat3 ipos, cfloat3 ivel, int cell, bufList buf)
//{
//	float dsq, c;
//	register float r2 = simData.r2;
//
//	cfloat3 dist, jvel;
//	float cmterm;
//
//	if (buf.mgridcnt[cell] == 0) return;
//
//	int cfirst = buf.mgridoff[cell];
//	int clast = cfirst + buf.mgridcnt[cell];
//
//	for (int j = cfirst; j < clast; j++) {
//		
//		jvel = buf.calcBuffer[j].veleval;
//
//		dist = (ipos - buf.displayBuffer[j].pos);
//		dist *= paramCarrier.simscale;
//		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
//
//		if (dsq < r2 && dsq > 0) {
//			dsq = sqrt(dsq);
//			c = (simData.psmoothradius - dsq);
//			cmterm = simData.spikykern * c * c / dsq * buf.calcBuffer[j].mass * buf.calcBuffer[j].dens;
//			//cmterm = simData.spikykern * c * c / dsq;
//			jvel = (jvel-ivel) * cmterm;
//			result[0] += jvel.x * dist.x;	result[1] += jvel.x * dist.y;	result[2] += jvel.x * dist.z;
//			result[3] += jvel.y * dist.x;	result[4] += jvel.y * dist.y;	result[5] += jvel.y * dist.z;
//			result[6] += jvel.z * dist.x;	result[7] += jvel.z * dist.y;	result[8] += jvel.z * dist.z;
//		}
//	}
//}


//__device__ void Drucker_Prager_Tensor(bufList buf, float* plastic, float* strain, float* tensor, int pid) {
//
//	// deviatoric stress
//	float lambda;
//
//	float s[9];
//	for (int i=0; i<9; i++) {
//		plastic[i] = 0; //initialize
//		s[i] = tensor[i];
//	}
//	float tr = (s[0]+s[4]+s[8])/3.0f;
//	s[0] -= tr;
//	s[4] -= tr;
//	s[8] -= tr;
//
//	//second invariant
//	float j2 = 0;
//	for (int i=0; i<9; i++)
//		j2 += 0.5*s[i]*s[i];
//	if (j2<=0.0f)
//		j2 = 0.0f;
//	else
//		j2 = sqrt(j2);
//
//	//first invariant
//	float i1 = tensor[0]+tensor[4]+tensor[8];
//
//	float a_phi = paramCarrier.a_phi;
//	float a_psi = paramCarrier.a_psi;
//	float k_c   = paramCarrier.k_c;
//
//	float f = j2 + a_phi*i1 - k_c;
//	if (f<=0.0f) { //no yield
//		return;
//	}
//	else {
//		//buf.mclr[pid] = COLORA(0, 0, 0, 1);
//		float strain_tr = strain[0]+strain[4]+strain[8];
//		float s_epsilon=0;
//		for (int i=0; i<9; i++) {
//			s_epsilon += s[i]*strain[i];
//		}
//
//		if (j2<0.0001f) {
//			lambda = 0.0f;
//			j2 = 1.0f;
//		}
//		else {
//			lambda = 3*a_phi*simData.solid_coK*strain_tr + simData.solid_coG/j2*s_epsilon;
//			lambda /= (9*a_phi*a_psi*simData.solid_coK + simData.solid_coG);
//		}
//	}
//	int id;
//	for (int ii=0; ii<3; ii++) {
//		for (int ij=0; ij<3; ij++) {
//			id = ii*3+ij;
//			plastic[id] = simData.solid_coG/j2*s[id] * lambda;
//			if (ii==ij)
//				plastic[id] += 3*a_psi*simData.solid_coK * lambda;
//		}
//	}
//
//
//}

//__global__ void ComputeSolidTensor_CUDA(bufList buf, int pnum) {
//	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
//	if (i >= pnum) return;
//	// Get search cell
//	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
//	uint gc = buf.mgcell[i];
//	if (gc == GRID_UNDEF) return;						// particle out-of-range
//	gc -= nadj;
//
//	cfloat3 ipos = buf.displayBuffer[i].pos;
//	cfloat3 ivel = buf.calcBuffer[i].vel;
//	float vgrad[9];
//	for (int ii=0; ii<9; ii++)
//		vgrad[ii]=0;
//
//	//Velocity Gradient
//	for (int c=0; c < simData.gridAdjCnt; c++) {
//		contributeVGrad(vgrad, i, ipos, ivel, gc+simData.gridAdj[c], buf);
//	}
//	//Strain
//	float epsilon[9], omega[9];
//	for (int ii=0; ii<3; ii++) {
//		for (int ij=0; ij<3; ij++) {
//			epsilon[ii*3+ij] = 0.5 * (vgrad[ii*3+ij] + vgrad[ij*3+ii]);
//			omega[ii*3+ij] = 0.5 * (vgrad[ii*3+ij] - vgrad[ij*3+ii]);
//		}
//	}
//	float strain = epsilon[0] + epsilon[4] + epsilon[8];
//
//	//Elastic Stress Update
//	float* stress = buf.MFtensor + i*9;
//	float* dev_stress = buf.MFtemptensor + i*9;
//
//	float d_sigma[9];
//	int id;
//	for (int ii=0; ii<3; ii++) {
//		for (int ij=0; ij<3; ij++) {
//			id = ii*3+ij;
//			d_sigma[id] = 0;
//			for (int ik=0; ik<3; ik++) {
//				d_sigma[id] += dev_stress[ii*3+ik]*omega[ij*3+ik];
//				d_sigma[id] += dev_stress[ik*3+ij]*omega[ii*3+ik];
//			}
//			if (ii==ij) {
//				d_sigma[id] += 2*simData.solid_coG * (epsilon[id] - strain/3);
//				//d_sigma[id] += simData.solid_coK * (strain);
//			}
//			else
//				d_sigma[id] += 2*simData.solid_coG * epsilon[id];
//		}
//	}
//
//	//Drucker-Prager Plastic Stress
//	//float plastic[9];
//
//	for (int ii=0; ii<9; ii++) {
//		dev_stress[ii] += d_sigma[ii] * simData.mf_dt;
//		stress[ii] = dev_stress[ii];
//	}
//	stress[0] -= buf.mpress[i];
//	stress[4] -= buf.mpress[i];
//	stress[8] -= buf.mpress[i];
//
//	//Drucker_Prager_Tensor(buf, plastic, epsilon, stress, i);
//
//	//for (int mi=0; mi<9; mi++) {
//	//	stress[mi] -= plastic[mi] * simData.mf_dt;
//	//}
//
//	//R-tensor
//	float* rtensor = buf.MFRtensor+i*9;
//	float sigma[3];
//	float u[9], v[9], svalue[9];
//
//	for (int ii=0; ii<9; ii++) {
//		u[ii] = stress[ii];
//		rtensor[ii] = 0;
//		svalue[ii] = 0;
//	}
//	svdecomp3(sigma, u, v, 0.000001f);
//	float res[9];
//
//	bool err = false;
//	for (int ii=0; ii<3; ii++) {
//		if (u[ii]*v[ii]<0) {
//			sigma[ii] *= -1;
//			v[ii]*=-1;
//			v[ii+3]*=-1;
//			v[ii+6]*=-1;
//		}
//		if (sigma[ii]<0)
//			sigma[ii] = 0;
//	}
//
//	float r[9];
//	transpose3(v, v);
//	svalue[0] = sigma[0];
//	svalue[4] = sigma[1];
//	svalue[8] = sigma[2];
//
//	multiply_matrix3(u, svalue, r);
//	multiply_matrix3(r, v, rtensor);
//
//	//if(i%100==0){
//	//print9("stress",i,stress);
//	//print9("r",i,r);
//	//print9("u",i,u);
//	//print9("v",i,v);
//	//print9("rtensor",i, rtensor);
//	//}
//
//	for (int ii=0; ii<9; ii++) {
//		rtensor[ii] = rtensor[ii] * paramCarrier.SFlipControl;
//	}
//
//}
//
//__device__ cfloat3 contributeSolidForce(cfloat3* force, int i, cfloat3 ipos, float* itensor, int cell, bufList buf)
//{
//	if (buf.mgridcnt[cell] == 0)
//		return make_cfloat3(0, 0, 0);
//	float dsq, c;
//	register float d2 = simData.psimscale * simData.psimscale;
//	register float r2 = simData.r2/d2;
//
//	cfloat3 dist, origin_dist;
//	float* jtensor;
//	float cmterm;
//	//	float massj,massi;
//	int j;
//	cfloat3 vij, xij;
//
//	cfloat3 result = make_cfloat3(0, 0, 0);
//
//	int cfirst = buf.mgridoff[cell];
//	int clast = cfirst + buf.mgridcnt[cell];
//
//	//massi = buf.mf_restmass[i];
//	float tensorij[9];
//	for (int k=0; k<9; k++)
//		tensorij[k] = 0;
//
//	int mftype = buf.MFtype[i];
//	int mftypej;
//	float rij[9];
//
//	for (int cndx = cfirst; cndx < clast; cndx++) {
//		j = buf.mgrid[cndx];
//		mftypej = buf.MFtype[j];
//
//		dist = (ipos - buf.mpos[j]);		// dist in cm
//		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
//		dist *= simData.psimscale;
//		origin_dist = -dist;
//		if (!(dsq < r2 && dsq > 0))
//			continue;
//
//		vij =  buf.mveleval[i] - buf.mveleval[j];
//		xij =  (ipos - buf.mpos[j])*simData.psimscale;
//
//		jtensor = buf.MFtensor + j*9;
//		//sigma force related to pressure
//
//		for (int k=0; k<9; k++) {
//			tensorij[k] = (buf.MFtensor[i*9+k]*buf.mdensity[i]*buf.mdensity[i] + buf.MFtensor[j*9+k]*buf.mdensity[j]*buf.mdensity[j]);
//			//tensorij[k] = (buf.MFtensor[i*9+k]+buf.MFtensor[j*9+k])*buf.mdensity[i]*buf.mdensity[j];
//		}
//		dsq = sqrt(dsq * d2);
//		c = (simData.psmoothradius - dsq);
//		cmterm = simData.spikykern * c * c / dsq * buf.mf_restmass[j];
//		dist *= cmterm;
//
//
//
//		//artificial viscosity
//		float xvprod = xij.x * vij.x + xij.y * vij.y + xij.z * vij.z;
//		float PIij;
//		if (xvprod < 0) {
//			float densityij = (1/buf.mdensity[i]+1/buf.mdensity[j])*0.5;
//			float ci = sqrt(4*simData.solid_coG/3*buf.mdensity[i] + simData.solid_coK*buf.mdensity[i]);
//			float cj = sqrt(4*simData.solid_coG/3*buf.mdensity[j] + simData.solid_coK*buf.mdensity[j]);
//			ci = (ci+cj)*0.5;
//			float phiij = simData.psmoothradius * xvprod/
//				((xij.x*xij.x + xij.y*xij.y + xij.z*xij.z) + 0.01* simData.r2);
//
//			PIij = (-simData.solid_coA *ci* phiij + simData.solid_coB * phiij*phiij)/densityij;
//
//			result.x -= PIij * dist.x;
//			result.y -= PIij * dist.y;
//			result.z -= PIij * dist.z;
//		}
//
//		//tensile instability
//		float Fij = simData.spikykernel * pow((simData.psmoothradius - dsq), 3)/paramCarrier.w_deltax;
//		Fij = pow(Fij, paramCarrier.fijn);
//		for (int ii=0; ii<9; ii++) {
//			rij[ii] = (buf.MFRtensor[i*9+ii]*buf.mdensity[i]*buf.mdensity[i] + buf.MFRtensor[j*9+ii]*buf.mdensity[j]*buf.mdensity[j])*Fij;
//			tensorij[ii] += rij[ii];
//		}
//
//		result.x += tensorij[0]*dist.x + tensorij[1]*dist.y + tensorij[2]*dist.z;
//		result.y += tensorij[3]*dist.x + tensorij[4]*dist.y + tensorij[5]*dist.z;
//		result.z += tensorij[6]*dist.x + tensorij[7]*dist.y + tensorij[8]*dist.z;
//
//		//if(i==0 && Fij>0.1){
//		//	//printf("h-dsq %f %f\n",simData.psmoothradius,dsq);
//		//	printf("%d fij %f %f\n",j, simData.spikykernel * pow((simData.psmoothradius - dsq), 3),dsq/simData.psimscale);
//		//	printf("wdeltax %f\n", paramCarrier.w_deltax);
//		//	print9("tensorij",j,tensorij);
//		//	printf("fij %f\n",Fij);
//		//	print9("rij",j,rij);
//		//}
//
//	}
//	*force += result;
//}
//
//
//
//__global__ void ComputeSolidForce_CUDA(bufList buf, int pnum) {
//	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
//	if (i >= pnum) return;
//	// Get search cell
//	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
//	uint gc = buf.mgcell[i];
//	if (gc == GRID_UNDEF) return;						// particle out-of-range
//	gc -= nadj;
//
//	cfloat3 ipos = buf.mpos[i];
//	cfloat3 ivel = buf.mveleval[i];
//	float* itensor = buf.MFtensor + i*9;
//	cfloat3 force = make_cfloat3(0, 0, 0);
//	//Force
//	for (int c=0; c < simData.gridAdjCnt; c++) {
//		contributeSolidForce(&force, i, ipos, itensor, gc+simData.gridAdj[c], buf);
//	}
//	buf.mforce[i] = force;
//}
//
//__device__ void contributeSolidDensity(uint i, float& res, uint cell, bufList buf) {
//	cfloat3 dist;
//	cfloat3 p = buf.mpos[i];
//	float dsq, c;
//	float massj;
//	register float d2 = simData.psimscale * simData.psimscale;
//	register float r2 = simData.r2/d2;
//
//	int j;
//
//	if (buf.mgridcnt[cell] == 0)
//		return;
//
//	int cfirst = buf.mgridoff[cell];
//	int clast = cfirst + buf.mgridcnt[cell];
//
//	for (int cndx = cfirst; cndx < clast; cndx++) {
//		j = buf.mgrid[cndx];
//		if (buf.misbound[j]==0)
//			massj = buf.mf_restmass[j];
//
//		dist = p - buf.mpos[j];
//		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
//		if (dsq < r2 && dsq > 0.0) {
//			c = (r2 - dsq)*d2;
//			res += c * c * c * massj;
//		}
//	}
//	return;
//}
//
//__global__ void ComputeDensity_CUDA(bufList buf, int pnum) {
//	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
//	if (i >= pnum) return;
//	// Get search cell
//	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
//	uint gc = buf.mgcell[i];
//	if (gc == GRID_UNDEF) return;						// particle out-of-range
//	gc -= nadj;
//
//	float sum=0;
//
//	for (int c=0; c < simData.gridAdjCnt; c++) {
//		contributeSolidDensity(i, sum, gc+simData.gridAdj[c], buf);
//	}
//
//	sum += simData.r2 * simData.r2 * simData.r2 * buf.mf_restmass[i];
//	sum = sum * simData.poly6kern;
//	buf.mdensity[i] = 1.0f/sum;
//	float dens = buf.mf_restdensity[i];
//
//	buf.mpress[i] = 2.5 * dens * (pow(sum / dens, 7.0f) - 1);
//}


/*----------------------------------

	accumulative for plastic

------------------------------------*/

__device__ void contributeVGrad(cmat3& result, int i, bufList buf, int cell)
{
	if (buf.mgridcnt[cell] == 0)
		return;

	float dsq, c;
	register float r2 = simData.r2;
	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 ivel = buf.calcBuffer[i].vel;

	cfloat3 dist, jvel;
	float cmterm;

	

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	for (int j = cfirst; j < clast; j++) {

		jvel = buf.calcBuffer[j].vel;

		dist = (ipos - buf.displayBuffer[j].pos);
		dist *= paramCarrier.simscale;
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);

		if (dsq < r2 && dsq > 0) {
			dsq = sqrt(dsq);
			c = (simData.psmoothradius - dsq);
			cmterm = simData.spikykern * c * c / dsq * buf.calcBuffer[j].mass * buf.calcBuffer[j].dens;
			//cmterm = simData.spikykern * c * c / dsq;
			jvel = (jvel - ivel) * cmterm;
			result[0][0] += jvel.x * dist.x;	result[0][1] += jvel.x * dist.y;	result[0][2] += jvel.x * dist.z;
			result[1][0] += jvel.y * dist.x;	result[1][1] += jvel.y * dist.y;	result[1][2] += jvel.y * dist.z;
			result[2][0] += jvel.z * dist.x;	result[2][1] += jvel.z * dist.y;	result[2][2] += jvel.z * dist.z;
		}
	}
}

__global__ void ComputeSolidTensor(bufList buf, int pnum) {

	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if(i>=pnum)
		return;
	if(buf.displayBuffer[i].type != TYPE_ELASTIC)
		return;

	//velocity gradient nabla v
	uint gc = buf.mgcell[i];
	if(gc==GRID_UNDEF)
		return;
	cmat3 velGrad;
	velGrad.Set(0.0f);
	
	for(int c=0; c<simData.gridAdjCnt; c++)
		contributeVGrad(velGrad, i, buf, gc+simData.gridAdj[c]);

	//deformation gradient F
	cmat3 tmp;
	mat3prod(velGrad, buf.calcBuffer[i].deformGrad, tmp);
	for (int k=0; k<9; k++) {
		buf.calcBuffer[i].deformGrad.data[k] += tmp.data[k] * paramCarrier.dt;
	}	
	tmp = buf.calcBuffer[i].deformGrad;

	//get strain epsilon
	cmat3 tmp2;
	mat3transpose(tmp, tmp2);
	cmat3 tmp3;
	mat3prod(tmp,tmp2,tmp3);
	tmp3[0][0] -= 1;
	tmp3[1][1] -= 1;
	tmp3[2][2] -= 1;

	//get Cauchy stress sigma
	for(int k=0; k<9; k++)
		tmp3.data[k] *= paramCarrier.solidK;
	buf.intmBuffer[i].stress = tmp3;

}

__device__ void ContributeSolidForce(cfloat3& result, int i, bufList buf, int cell) {
	if (buf.mgridcnt[cell] == 0)
		return;

	float dsq, c;
	register float r2 = simData.r2;
	cfloat3 ipos = buf.displayBuffer[i].pos;

	cfloat3 dist;
	float cmterm;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	cmat3 stressij;
	float idens,jdens;
	idens = buf.calcBuffer[i].dens;

	for (int j = cfirst; j < clast; j++) {

		dist = (ipos - buf.displayBuffer[j].pos);
		dist *= paramCarrier.simscale;
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);

		if (dsq < r2 && dsq > 0) {
			dsq = sqrt(dsq);
			c = (simData.psmoothradius - dsq);
			jdens = buf.calcBuffer[j].dens;

			if (buf.displayBuffer[j].type == TYPE_ELASTIC){
				
				cmterm = simData.spikykern * c * c / dsq * buf.calcBuffer[j].mass * jdens * idens;

				for(int k=0; k<9; k++)
					stressij.data[k] = (buf.intmBuffer[i].stress.data[k] + buf.intmBuffer[j].stress.data[k]);

				cfloat3 tmp;
				mvprod(stressij, dist, tmp);
				tmp *= cmterm;
				result += tmp;

				/*if (buf.calcBuffer[i].bornid==0) {
					printf("stressij\n");
					stressij.Print();
					printf("\n");
				}*/
			}
		}
	}
}

__global__ void ComputeSolidForce(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i>=pnum)
		return;
	if (buf.displayBuffer[i].type != TYPE_ELASTIC)
		return;

	uint gc = buf.mgcell[i];
	if (gc==GRID_UNDEF)
		return;
	
	cfloat3 accel(0,0,0);
	for (int c=0; c<simData.gridAdjCnt; c++)
		ContributeSolidForce(accel, i, buf, gc+simData.gridAdj[c]);
	
	buf.calcBuffer[i].accel = accel;

	/*if (buf.calcBuffer[i].bornid==0) {
		printf("%f %f %f\n", accel.x, accel.y, accel.z);
		printf("\n");
	}*/
}


/*----------------------------------

	hyperelastic material

-----------------------------------*/


__device__ void contributeDeformGrad(cmat3& res, int i, bufList buf, int cell) {

}

//with reference shape
__global__ void ComputeSolidTensor_X(bufList buf, int pnum) {

}

__device__ void contributeSolidForce_X(bufList buf, int pnum) {

}

__global__ void  ComputeSolidForce_X(bufList buf, int pnum) {

}

__global__ void ComputeInvA(bufList buf, int pnum) {
	


}