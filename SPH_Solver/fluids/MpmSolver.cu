
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>


#include "fluid_system_kern.cuh"
#include "thrust\device_vector.h"	//thrust libs
#include "thrust\sort.h"
#include "fluidMath.cuh"

//extern __device__ FluidParams	simData;
extern __device__ ParamCarrier paramCarrier;

//Begin MPM Solution

__device__ void TensorProd(float* a, float* b, float* c){
    //a vector
    //b vector
    //c tensor
    for (int i=0; i<3; i++)
    for (int j=0; j<3; j++){
        c[i*3+j] = a[i] * b[j];
    }
}

__inline__ __device__ void TensorProd(cfloat3& a, cfloat3& b, cmat3& c) {
	c[0][0] = a.x*b.x; c[0][1] = a.x*b.y; c[0][2] = a.x*b.z;
	c[1][0] = a.y*b.x; c[1][1] = a.y*b.y; c[1][2] = a.y*b.z;
	c[2][0] = a.z*b.x; c[2][1] = a.z*b.y; c[2][2] = a.z*b.z;
}

__device__ int GetMpmIdx(int xid, int yid, int zid) {
	cint3& res = paramCarrier.mpmRes;
	return yid*res.x*res.z + zid*res.x + xid;
}

__global__ void initMpm(bufList buf, int mpmSize){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;
    if ( i >= mpmSize )
        return;

	cint3 mpmres = paramCarrier.mpmRes;
    int yid = i / mpmres.x / mpmres.z;
    int xid = i % mpmres.x;
    int zid = i / mpmres.x % mpmres.z;

    cint3 gridRes = paramCarrier.gridres;
    cint3 gridScan = gridRes;
	gridScan.x -= 1;
	gridScan.y -= 1;
	gridScan.z -= 1;


    cfloat3 pos;
	cfloat3 dx(xid,yid,zid);
	pos = paramCarrier.mpmXmin + dx * paramCarrier.mpmcellsize;
    buf.mpmPos[i] = pos;


    cfloat3 gcf = (pos - paramCarrier.gridmin) * paramCarrier.gridIdfac;
    cint3 gc = cint3 ( int ( gcf.x ), int ( gcf.y ), int ( gcf.z ) );    
    uint gs = (gc.y * gridRes.z + gc.z)*gridRes.x + gc.x;
    if (gc.x >= 1 && gc.x <= gridScan.x && gc.y >= 1 && gc.y <= gridScan.y && gc.z >= 1 && gc.z <= gridScan.z) {
		buf.mpmGid[i] = gs;
    }
    else{
        buf.mpmGid[i] = GRID_UNDEF;
		printf("%f %f %f %d %d %d\n", pos.x, pos.y, pos.z, xid, yid, zid);  
    }

}


__global__ void MpmColorTest(bufList buf, int mpmSize) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= mpmSize)
		return;

	if (buf.mpmPos[i].x<0) {
		uint gc = buf.mpmGid[i];
		if (gc==GRID_UNDEF) {
			printf("mpm nodex index error.\n");
			return;
		}

		for (int c=0; c<paramCarrier.neighbornum; c++) {
			int cell = gc + paramCarrier.neighborid[c];
			
			if(buf.mgridcnt[cell]==0)
				continue;

			int cfirst = buf.mgridoff[cell];
			int clast = cfirst + buf.mgridcnt[cell];

			for (int j=cfirst; j<clast; j++) {
				buf.displayBuffer[j].color.Set(1,1,1,1);
				
				//if(i%10==0)
					
			}
		}
	}
}

__global__ void MpmGetMomentum(bufList buf, int mpmSize) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= mpmSize)
		return;

	
	uint gc = buf.mpmGid[i];
	cfloat3 pos = buf.mpmPos[i];
	
	cfloat3 posj; 
	cfloat3 posij;
	float dist;
	float weightij;

	cfloat3 velsum(0,0,0);
	float weightsum=0;
	float cellsize = paramCarrier.mpmcellsize;

	for (int c=0; c<paramCarrier.neighbornum; c++) {
		int cell = gc + paramCarrier.neighborid[c];

		if (buf.mgridcnt[cell]==0)
			continue;

		int cfirst = buf.mgridoff[cell];
		int clast = cfirst + buf.mgridcnt[cell];

		for (int j=cfirst; j<clast; j++) {
			posj = buf.displayBuffer[j].pos; //particle
			posij = pos - posj;
			
			if(abs(posij.x)>=cellsize || abs(posij.y)>=cellsize || abs(posij.z)>=cellsize)
				continue;

			//trilinear interpolation
			weightij = (1-abs(posij.x)/cellsize) * (1-abs(posij.y)/cellsize) * (1-abs(posij.z)/cellsize);
			velsum += buf.calcBuffer[j].vel * weightij;
			weightsum += weightij;
		}
	}
	if(weightsum > EPSILON){
		buf.mpmVel[i] = velsum / weightsum;
		//printf("%f %f %f\n",buf.mpmVel[i].x, buf.mpmVel[i].y, buf.mpmVel[i].z);
	}
	else
		buf.mpmVel[i].Set(0,0,0);
}

__inline__ __device__ float Spline1D(float x) {
	x = abs(x);
	if(x>2)
		return 0;
	else if (x>1) {
		return -0.166666f*x*x*x + x*x - 2*x + 1.33333f;
	}
	else {
		return 0.5f*x*x*x - x*x + 0.666666f;
	}
}



__inline__ __device__ float SplineWeight(cfloat3 xij) {
	return Spline1D(xij.x) * Spline1D(xij.y) * Spline1D(xij.z);
}

__inline__ __device__ float NablaSpline1D(float x) {
	if(x<=-2)
		return 0;
	else if(x<=-1)
		return 0.5f*x*x + 2*x + 2;
	else if(x<=0)
		return -1.5f*x*x - 2*x;
	else if(x<1)
		return 1.5f*x*x - 2*x;
	else if(x<2)
		return -0.5f*x*x + 2*x - 2;
	else
		return 0;
}

__inline__ __device__ cfloat3 NablaSpline(cfloat3 xij) {
	cfloat3 Nx(Spline1D(xij.x), Spline1D(xij.y), Spline1D(xij.z));
	cfloat3 res;
	res.x = NablaSpline1D(xij.x)*Nx.y*Nx.z;
	res.y = NablaSpline1D(xij.y)*Nx.x*Nx.z;
	res.z = NablaSpline1D(xij.z)*Nx.x*Nx.y;
	return res;
}

__global__ void MpmParticleToGrid(bufList buf, int mpmSize) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= mpmSize)
		return;

	uint gc = buf.mpmGid[i];
	cfloat3 pos = buf.mpmPos[i];

	cfloat3 posj;
	cfloat3 posij;
	float dist;
	
	cfloat3 momentumSum(0, 0, 0);
	float massSum=0;
	float weight;
	float cellsize = paramCarrier.mpmcellsize;

	for (int c=0; c<paramCarrier.neighbornum; c++) {
		int cell = gc + paramCarrier.neighborid[c];

		if (buf.mgridcnt[cell]==0)
			continue;

		int cfirst = buf.mgridoff[cell];
		int clast = cfirst + buf.mgridcnt[cell];

		for (int j=cfirst; j<clast; j++) {
			if(buf.displayBuffer[j].type==TYPE_BOUNDARY)
				continue;

			posj = buf.displayBuffer[j].pos; //particle
			posij = (pos - posj)/cellsize;

			weight = SplineWeight(posij);
			massSum += weight * buf.calcBuffer[j].mass;
			momentumSum += buf.calcBuffer[j].vel * buf.calcBuffer[j].mass * weight;
		}
	}
	if (massSum > 0.000000001f) {
		buf.mpmVel[i] = momentumSum / massSum;
		buf.mpmMass[i] =  massSum;
	}
	else {
		buf.mpmVel[i].Set(0, 0, 0);
		buf.mpmMass[i] = 0;
	}
}

__device__ void Drucker_Prager_Tensor(bufList buf, float* plastic, float* strain, float* tensor, int pid) {

	// deviatoric stress
	float lambda;

	float s[9];
	for (int i=0; i<9; i++) {
		plastic[i] = 0; //initialize
		s[i] = tensor[i];
	}
	float tr = (s[0]+s[4]+s[8])/3.0f;
	s[0] -= tr;
	s[4] -= tr;
	s[8] -= tr;

	//second invariant
	float j2 = 0;
	for (int i=0; i<9; i++)
		j2 += 0.5*s[i]*s[i];
	if (j2<=0.0f)
		j2 = 0.0f;
	else
		j2 = sqrt(j2);

	//first invariant
	float i1 = tensor[0]+tensor[4]+tensor[8];

	float a_phi = paramCarrier.a_phi;
	float a_psi = paramCarrier.a_psi;
	float k_c   = paramCarrier.k_c;

	float f = j2 + a_phi*i1 - k_c;
	if (f<=0.0f) { //no yield
		return;
	}
	else {
		//buf.mclr[pid] = COLORA(0, 0, 0, 1);
		float strain_tr = strain[0]+strain[4]+strain[8];
		float s_epsilon=0;
		for (int i=0; i<9; i++) {
			s_epsilon += s[i]*strain[i];
		}

		if (j2<0.0001f) {
			lambda = 0.0f;
			j2 = 1.0f;
		}
		else {
			lambda = 3*a_phi*paramCarrier.solidK*strain_tr + paramCarrier.solidG/j2*s_epsilon;
			lambda /= (9*a_phi*a_psi*paramCarrier.solidK + paramCarrier.solidG);
		}
	}
	int id;
	for (int ii=0; ii<3; ii++) {
		for (int ij=0; ij<3; ij++) {
			id = ii*3+ij;
			plastic[id] = paramCarrier.solidG/j2*s[id] * lambda;
			if (ii==ij)
				plastic[id] += 3*a_psi*paramCarrier.solidK * lambda;
		}
	}


}

__device__ cmat3& Drucker_Prager_Tensor(bufList buf, cmat3& epsilon, cmat3& sigma, int pid) {

	// deviatoric stress
	float lambda;
	cmat3 s, plastic;
	plastic.Set(0.0f);
	s = sigma;
	
	float trsigma = (s[0][0]+s[1][1]+s[2][2])/3.0f;
	s[0][0] -= trsigma;
	s[1][1] -= trsigma;
	s[2][2] -= trsigma;

	//second invariant
	float j2 = 0;
	for (int i=0; i<9; i++)
		j2 += 0.5*s.data[i]*s.data[i];
	
	if (j2 <= 0.0000001f) {
		return plastic;
	}
	else
		j2 = sqrt(j2);

	//first invariant
	float i1 = sigma[0][0] + sigma[1][1] + sigma[2][2];

	float a_phi = paramCarrier.a_phi;
	float a_psi = paramCarrier.a_psi;
	float k_c   = paramCarrier.k_c;

	float f = j2 + a_phi*i1 - k_c;
	if (f<=0.0f) { //no yield
		//buf.displayBuffer[pid].color.Set(0,1,0,1);
		return plastic;
	}
	else {
		float strain_tr = epsilon[0][0] + epsilon[1][1] + epsilon[2][2];
		float s_epsilon=0;
		for (int i=0; i<9; i++) {
			s_epsilon += s.data[i]*epsilon.data[i];
		}

		lambda = 3*a_phi*paramCarrier.solidK*strain_tr + paramCarrier.solidG/j2*s_epsilon;
		lambda /= (9*a_phi*a_psi*paramCarrier.solidK + paramCarrier.solidG);
		//buf.displayBuffer[pid].color.Set(1,0,0,1);
	}
	
	for (int ii=0; ii<3; ii++) {
		for (int ij=0; ij<3; ij++) {
			
			plastic[ii][ij] = paramCarrier.solidG/j2 * s[ii][ij] * lambda;
			if (ii==ij)
				plastic[ii][ij] += 3*a_psi*paramCarrier.solidK * lambda;
		}
	}


}

__global__ void MpmParticleStressFiniteStrain(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum)
		return;
	if(buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;
	//particle vel gradient
	
	cfloat3 pos = buf.displayBuffer[i].pos;
	cfloat3 xj;
	cfloat3 xij;
	float dist;

	cfloat3 velSum(0, 0, 0);
	float weight;
	float cellsize = paramCarrier.mpmcellsize;
	float realcellsize = cellsize * paramCarrier.simscale;
	//get particle mpmgrid cellid
	int xid, yid, zid;
	xid = (pos.x - paramCarrier.mpmXmin.x)/cellsize;
	yid = (pos.y - paramCarrier.mpmXmin.y)/cellsize;
	zid = (pos.z - paramCarrier.mpmXmin.z)/cellsize;
	int mpmid = GetMpmIdx(xid, yid, zid);
	cint3& res = paramCarrier.mpmRes;
	cfloat3 nablaWij;
	cmat3 v_x_nw;
	cmat3 velGrad;
	velGrad.Set(0.0f);
	float densityp = 0.0f;

	int tx, ty, tz, tid;
	for (int ii=-1; ii<3; ii++) {
		for (int ij=-1; ij<3; ij++) {
			for (int ik=-1; ik<3; ik++) {
				//4*4*4
				tx = xid+ii;
				ty = yid+ij;
				tz = zid+ik;
				if (tx<0 || tx>=res.x || ty<0 || ty>=res.y || tz<0 || tz>=res.z)
					continue;
				tid = GetMpmIdx(tx, ty, tz);
				xj = buf.mpmPos[tid];
				xij = (pos - xj)/cellsize;
				
				weight = SplineWeight(xij);
				//velSum += buf.mpmVel[tid]*weight;

				nablaWij = NablaSpline(xij);
				TensorProd(buf.mpmVel[tid], nablaWij, v_x_nw);
				mat3add(velGrad, v_x_nw, velGrad);

				densityp += buf.mpmMass[tid] * weight / realcellsize / realcellsize / realcellsize;
			}
		}
	}
	buf.calcBuffer[i].dens = densityp;

	//stress update
	//deformation gradient F
	cmat3 tmp;
	mat3prod(velGrad, buf.calcBuffer[i].deformGrad, tmp);
	for (int k=0; k<9; k++) {
		buf.calcBuffer[i].deformGrad.data[k] += tmp.data[k] * paramCarrier.dt;
	}
	tmp = buf.calcBuffer[i].deformGrad;

	//get strain epsilon //epsilon = J J^T - I
	cmat3 tmp2;
	mat3transpose(tmp, tmp2);
	cmat3 tmp3;
	mat3prod(tmp2, tmp, tmp3);
	tmp3[0][0] -= 1;
	tmp3[1][1] -= 1;
	tmp3[2][2] -= 1;

	//get Cauchy stress sigma
	for (int k=0; k<9; k++)
		tmp3.data[k] *= paramCarrier.solidK;
	buf.stress[i] = tmp3;

	//printf("stress %f %f %f\n %f %f %f\n %f %f %f\n", tmp3[0][0], tmp3[0][1], tmp3[0][2], 
	//	tmp3[1][0], tmp3[1][1], tmp3[1][2],
	//	tmp3[2][0], tmp3[2][1], tmp3[2][2]);
}

__global__ void MpmParticleStress(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum)
		return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;
	
	//particle vel gradient

	cfloat3 pos = buf.displayBuffer[i].pos;
	cfloat3 xj;
	cfloat3 xij;
	float dist;

	cfloat3 velSum(0, 0, 0);
	float weight;
	float cellsize = paramCarrier.mpmcellsize;
	float realcellsize = cellsize * paramCarrier.simscale;
	//get particle mpmgrid cellid
	int xid, yid, zid;
	xid = (pos.x - paramCarrier.mpmXmin.x)/cellsize;
	yid = (pos.y - paramCarrier.mpmXmin.y)/cellsize;
	zid = (pos.z - paramCarrier.mpmXmin.z)/cellsize;
	int mpmid = GetMpmIdx(xid, yid, zid);
	cint3& res = paramCarrier.mpmRes;
	cfloat3 nablaWij;
	cmat3 v_x_nw;
	cmat3 velGrad;
	velGrad.Set(0.0f);
	float densityp = 0.0f;

	int tx, ty, tz, tid;
	for (int ii=-1; ii<3; ii++) {
		for (int ij=-1; ij<3; ij++) {
			for (int ik=-1; ik<3; ik++) {
				//4*4*4
				tx = xid+ii;
				ty = yid+ij;
				tz = zid+ik;
				if (tx<0 || tx>=res.x || ty<0 || ty>=res.y || tz<0 || tz>=res.z)
					continue;
				tid = GetMpmIdx(tx, ty, tz);
				xj = buf.mpmPos[tid];
				xij = (pos - xj)/cellsize;

				weight = SplineWeight(xij);
			
				nablaWij = NablaSpline(xij);
				TensorProd(buf.mpmVel[tid], nablaWij, v_x_nw);
				mat3add(velGrad, v_x_nw, velGrad);

				densityp += buf.mpmMass[tid] * weight / realcellsize / realcellsize / realcellsize;
			}
		}
	}
	buf.calcBuffer[i].dens = densityp;

	//Cauchy Stress Tensor - infinitesimal strain theory
	cmat3 depsilon, omega;
	cmat3 vgradT;
	mat3transpose(velGrad, vgradT);
	for (int ii=0; ii<9; ii++) {
		depsilon.data[ii] = (velGrad.data[ii]+vgradT.data[ii])*0.5f;
		omega.data[ii] = (velGrad.data[ii]-vgradT.data[ii])*0.5f;
	}
	float trEpsilon = (depsilon[0][0] + depsilon[1][1] + depsilon[2][2])/3.0f;
	//depsilon[0][0] -= trEpsilon;
	//depsilon[1][1] -= trEpsilon;
	//depsilon[2][2] -= trEpsilon;
	cmat3 dsigma;
	cmat3& sigma = buf.calcBuffer[i].stress;

	for (int ii=0; ii<3; ii++) {
		for (int ij=0; ij<3; ij++) {
			dsigma[ii][ij] = 2*paramCarrier.solidG*depsilon[ii][ij];
			for (int ik=0; ik<3; ik++) {
				dsigma[ii][ij] += omega[ii][ik]*sigma[ik][ij] - sigma[ii][ik]*omega[ik][ij];
			}
			if(ii==ij)
				dsigma[ii][ii] += paramCarrier.solidK * trEpsilon * 3.0f - 2*paramCarrier.solidG*trEpsilon;
		}
	}
	
	cmat3 plastic = Drucker_Prager_Tensor(buf, depsilon, sigma, i);

	for (int ii=0; ii<9; ii++) {
		sigma.data[ii] += dsigma.data[ii] * paramCarrier.dt;
		sigma.data[ii] -= plastic.data[ii] * paramCarrier.dt;
	}

	buf.stress[i] = sigma;
}


__device__ void GridCollisionDetection(cfloat3& pos, cfloat3& vel) {
	// Get particle vars
	register cfloat3 accel, norm;
	register float diff, adj, speed;
	float time = 0;
	float ss = paramCarrier.simscale;

	accel = cfloat3(0, 0, 0);

	// Soft Boundaries
	// Y-axis
	diff = -(pos.y - paramCarrier.softminx.y) * ss;
	if (diff > EPSILON) {
		norm = cfloat3(0, 1.0, 0);
		if (dot(norm, vel)<0) {
			vel.Set(0,0,0);
		}
	}
	diff = -(pos.x - paramCarrier.softminx.x) * ss;
	if (diff > EPSILON) {
		norm = cfloat3(1.0, 0.0, 0);
		if (dot(norm, vel)<0) {
			vel.Set(0, 0, 0);
		}
	}
	diff = -(pos.z - paramCarrier.softminx.z) * ss;
	if (diff > EPSILON) {
		norm = cfloat3(0, 0.0, 1.0);
		if (dot(norm, vel)<0) {
			vel.Set(0, 0, 0);
		}
	}
	diff = (pos.x - paramCarrier.softmaxx.x) * ss;
	if (diff > EPSILON) {
		norm = cfloat3(-1.0, 0.0, 0);
		if (dot(norm, vel)<0) {
			vel.Set(0, 0, 0);
		}
	}
	diff = (pos.z - paramCarrier.softmaxx.z) * ss;
	if (diff > EPSILON) {
		norm = cfloat3(0, 0.0, -1.0);
		if (dot(norm, vel)<0) {
			vel.Set(0, 0, 0);
		}
	}


	//vel = accel*paramCarrier.dt + vel;				// v(t+1/2) = v(t-1/2) + a(t) dt
}

__global__ void MpmNodeUpdate(bufList buf, int mpmSize) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= mpmSize)
		return;
	if(buf.mpmMass[i]<0.000000001f)
		return;

	//node force
	uint gc = buf.mpmGid[i];
	cfloat3 pos = buf.mpmPos[i];

	cfloat3 posj;
	cfloat3 posij;
	float dist;

	cfloat3 forceSum(0, 0, 0);
	float cellsize = paramCarrier.mpmcellsize;
	cfloat3 nablaWij;
	cfloat3 sigma_nw;

	for (int c=0; c<paramCarrier.neighbornum; c++) {
		int cell = gc + paramCarrier.neighborid[c];

		if (buf.mgridcnt[cell]==0)
			continue;

		int cfirst = buf.mgridoff[cell];
		int clast = cfirst + buf.mgridcnt[cell];

		for (int j=cfirst; j<clast; j++) {
			if (buf.displayBuffer[j].type==TYPE_BOUNDARY)
				continue;

			posj = buf.displayBuffer[j].pos; //particle
			posij = (pos - posj)/cellsize;
			nablaWij = NablaSpline(posij);
			mvprod(buf.stress[j], nablaWij, sigma_nw);
			forceSum += sigma_nw / buf.calcBuffer[j].dens;
		}
	}
	forceSum += paramCarrier.gravity;
	//node velocity
	//printf("%f %f %f\n",forceSum.x, forceSum.y, forceSum.z);
	buf.mpmVel[i] += forceSum * paramCarrier.dt;
	GridCollisionDetection(buf.mpmPos[i], buf.mpmVel[i]);

}

__global__ void MpmParticleUpdate(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum)
		return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY) {
		return;
	}

	cfloat3 pos = buf.displayBuffer[i].pos;

	cfloat3 xj;
	cfloat3 xij;
	float dist;

	cfloat3 velSum(0, 0, 0);
	float weight;
	float cellsize = paramCarrier.mpmcellsize;

	//get particle mpmgrid cellid
	int xid,yid,zid;
	xid = (pos.x - paramCarrier.mpmXmin.x)/cellsize;
	yid = (pos.y - paramCarrier.mpmXmin.y)/cellsize;
	zid = (pos.z - paramCarrier.mpmXmin.z)/cellsize;
	int mpmid = GetMpmIdx(xid,yid,zid);
	cint3& res = paramCarrier.mpmRes;

	int tx,ty,tz,tid;
	for (int ii=-1; ii<3; ii++) {
		for (int ij=-1; ij<3; ij++) {
			for (int ik=-1; ik<3; ik++) {
				//4*4*4
				tx = xid+ii;
				ty = yid+ij;
				tz = zid+ik;
				if(tx<0 || tx>=res.x || ty<0 || ty>=res.y || tz<0 || tz>=res.z)
					continue;
				tid = GetMpmIdx(tx,ty,tz);
				xj = buf.mpmPos[tid];
				xij = (pos - xj)/cellsize;
				weight = SplineWeight(xij);
				velSum += buf.mpmVel[tid]*weight;
			}
		}
	}
	buf.calcBuffer[i].vel = velSum;
	/*if (buf.calcBuffer[i].bornid==0) {
		printf("pvel: %f %f %f\n",velSum.x, velSum.y, velSum.z);
	}*/
	//position integration
	buf.displayBuffer[i].pos += buf.calcBuffer[i].vel * paramCarrier.dt / paramCarrier.simscale;
}

__device__ void printmat(const char* str, cmat3& m) {
	printf("%s %.10f %.10f %.10f\n %.10f %.10f %.10f\n %.10f %.10f %.10f\n",str, m.data[0],m.data[1],m.data[2], m.data[3],m.data[4],m.data[5], m.data[6],m.data[7],m.data[8]);
}

__global__ void MpmParticleDp(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum)
		return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY) {
		return;
	}

	cfloat3 pos = buf.displayBuffer[i].pos;

	cfloat3 xj;
	cfloat3 xij;
	float dist;

	float weight;
	float cellsize = paramCarrier.mpmcellsize;
	cfloat3 xij_real;
	cmat3 xx;

	//get particle mpmgrid cellid
	int xid, yid, zid;
	xid = (pos.x - paramCarrier.mpmXmin.x)/cellsize;
	yid = (pos.y - paramCarrier.mpmXmin.y)/cellsize;
	zid = (pos.z - paramCarrier.mpmXmin.z)/cellsize;
	int mpmid = GetMpmIdx(xid, yid, zid);
	cint3& res = paramCarrier.mpmRes;
	cmat3 dp;
	dp.Set(0.0f);

	int tx, ty, tz, tid;
	for (int ii=-1; ii<3; ii++) {
		for (int ij=-1; ij<3; ij++) {
			for (int ik=-1; ik<3; ik++) {
				//4*4*4
				tx = xid+ii;
				ty = yid+ij;
				tz = zid+ik;
				if (tx<0 || tx>=res.x || ty<0 || ty>=res.y || tz<0 || tz>=res.z)
					continue;
				tid = GetMpmIdx(tx, ty, tz);
				xj = buf.mpmPos[tid];
				xij = (pos - xj)/cellsize;
				weight = SplineWeight(xij);
				
				xij_real = (pos - xj)*paramCarrier.simscale;
				TensorProd(xij_real, xij_real, xx);
				for (int iz=0; iz<9; iz++) {
					dp.data[iz] += xx.data[iz]*weight;
				}
			}
		}
	}

	//get inv Dp
	buf.invDp[i] = dp.Inv();
	if (buf.calcBuffer[i].bornid<=2) {
		printmat("dp",dp);
		printmat("inv dp", buf.invDp[i]);
	}
}

__global__ void MpmParticleToGrid_APIC(bufList buf, int mpmSize) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= mpmSize)
		return;

	//facility
	uint gc = buf.mpmGid[i];
	cfloat3 pos = buf.mpmPos[i];
	cfloat3 posj;
	cfloat3 posij;
	float dist;
	float weight;
	float cellsize = paramCarrier.mpmcellsize;

	//output
	cfloat3 velBD;
	cfloat3 momentumSum(0, 0, 0);
	float massSum=0;
	
	//constant
	float iDp;
	iDp = cellsize*cellsize * paramCarrier.simscale*paramCarrier.simscale / 3.0f;
	iDp = 1.0f / iDp;
	//if(i==0)
	//	printf("%f",iDp);
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		int cell = gc + paramCarrier.neighborid[c];

		if (buf.mgridcnt[cell]==0)
			continue;

		int cfirst = buf.mgridoff[cell];
		int clast = cfirst + buf.mgridcnt[cell];

		for (int j=cfirst; j<clast; j++) {
			if (buf.displayBuffer[j].type==TYPE_BOUNDARY)
				continue;

			posj = buf.displayBuffer[j].pos; //particle
			posij = (pos - posj)/cellsize;

			weight = SplineWeight(posij);
			massSum += weight * buf.calcBuffer[j].mass;
			
			
			cmat3 B = buf.calcBuffer[j].B;
			for (int iz=0; iz<9; iz++) {
				B.data[iz] *= iDp;
			}
			cfloat3 xij;
			xij = (pos - posj)*paramCarrier.simscale;
			mvprod(B, xij, velBD);
			velBD += buf.calcBuffer[j].vel;
			momentumSum += velBD * buf.calcBuffer[j].mass * weight;
		}
	}

	//check zero
	if (massSum > 0.000000001f) {
		buf.mpmVel[i] = momentumSum / massSum;
		buf.mpmMass[i] =  massSum;
	}
	else {
		buf.mpmVel[i].Set(0, 0, 0);
		buf.mpmMass[i] = 0;
	}
}

__global__ void MpmParticleUpdate_APIC(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum)
		return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY) {
		return;
	}

	//facility
	cfloat3 pos = buf.displayBuffer[i].pos;
	cfloat3 xj;
	cfloat3 xij;
	float dist;
	float weight;
	float cellsize = paramCarrier.mpmcellsize;

	//output
	cfloat3 velSum(0, 0, 0);
	cmat3& B = buf.calcBuffer[i].B;
	B.Set(0.0f);

	//get particle mpmgrid cellid
	int xid, yid, zid;
	xid = (pos.x - paramCarrier.mpmXmin.x)/cellsize;
	yid = (pos.y - paramCarrier.mpmXmin.y)/cellsize;
	zid = (pos.z - paramCarrier.mpmXmin.z)/cellsize;
	int mpmid = GetMpmIdx(xid, yid, zid);
	cint3& res = paramCarrier.mpmRes;

	int tx, ty, tz, tid;
	for (int ii=-1; ii<3; ii++) {
		for (int ij=-1; ij<3; ij++) {
			for (int ik=-1; ik<3; ik++) {
				//4*4*4
				tx = xid+ii;
				ty = yid+ij;
				tz = zid+ik;
				if (tx<0 || tx>=res.x || ty<0 || ty>=res.y || tz<0 || tz>=res.z)
					continue;
				tid = GetMpmIdx(tx, ty, tz);
				xj = buf.mpmPos[tid];
				xij = (pos - xj)/cellsize;
				weight = SplineWeight(xij);
				
				velSum += buf.mpmVel[tid]*weight;
				cmat3 vx;
				cfloat3 xij = (xj-pos)*paramCarrier.simscale;
				TensorProd(buf.mpmVel[tid], xij, vx);
				for (int iz=0; iz<9; iz++) {
					B.data[iz] += weight * vx.data[iz];
				}
			}
		}
	}
	buf.calcBuffer[i].vel = velSum;
	/*if (buf.calcBuffer[i].bornid==0) {
	printf("pvel: %f %f %f\n",velSum.x, velSum.y, velSum.z);
	}*/
	//position integration
	buf.displayBuffer[i].pos += buf.calcBuffer[i].vel * paramCarrier.dt / paramCarrier.simscale;
}


/*
__device__ void ShapeFunctionValue(cfloat3& xI, cfloat3& xp, float& val){
    int nxI, nyI, nzI;//natural coordinates
    cfloat3 center;
    if (abs(xp.x-xI.x)>simData.mpmSpacing ||
        abs(xp.y-xI.y)>simData.mpmSpacing ||
        abs(xp.z-xI.z)>simData.mpmSpacing){
        val = 0; 
        return;
    }

    if (xp.x > xI.x){
        nxI = -1;
        center.x = xI.x + simData.mpmSpacing * 0.5;
    }
    else{
        nxI = 1;
        center.x = xI.x - simData.mpmSpacing * 0.5;
    }
       
    if (xp.y  > xI.y){
        nyI = -1;
        center.y = xI.y + simData.mpmSpacing * 0.5;
    }
    else{
        nyI = 1;
        center.y = xI.y - simData.mpmSpacing * 0.5;
    }
    if (xp.z > xI.z){
        nzI = -1;
        center.z = xI.z + simData.mpmSpacing * 0.5;
    }
    else{
        nzI = 1;
        center.z = xI.z - simData.mpmSpacing * 0.5;
    }

    cfloat3 natCod = (xp-center)*2/simData.mpmSpacing;
    val = 0.125 * (1+natCod.x*nxI) * (1+natCod.y*nyI) * (1+natCod.z*nzI);
}


__device__ void ShapeFunctionGradValue(cfloat3& xI, cfloat3& xp, float* nGrad){
    int nxI, nyI, nzI;//natural coordinates
    cfloat3 center;
    if (abs(xp.x-xI.x)>simData.mpmSpacing ||
        abs(xp.y-xI.y)>simData.mpmSpacing ||
        abs(xp.z-xI.z)>simData.mpmSpacing){
        nGrad[0] = 0;
        nGrad[1] = 0;
        nGrad[2] = 0;
        return;
    }

    if (xp.x > xI.x){
        nxI = -1;
        center.x = xI.x + simData.mpmSpacing * 0.5;
    }
    else{
        nxI = 1;
        center.x = xI.x - simData.mpmSpacing * 0.5;
    }
       
    if (xp.y  > xI.y){
        nyI = -1;
        center.y = xI.y + simData.mpmSpacing * 0.5;
    }
    else{
        nyI = 1;
        center.y = xI.y - simData.mpmSpacing * 0.5;
    }
    if (xp.z > xI.z){
        nzI = -1;
        center.z = xI.z + simData.mpmSpacing * 0.5;
    }
    else{
        nzI = 1;
        center.z = xI.z - simData.mpmSpacing * 0.5;
    }

    cfloat3 natCod = (xp-center)*2/simData.mpmSpacing;
    float dxN[3];

    nGrad[0] = 0.25 * nxI * (1+nyI*natCod.y) * (1+nzI*natCod.z) / simData.mpmSpacing;
    nGrad[1] = 0.25 * nyI * (1+nxI*natCod.x) * (1+nzI*natCod.z) / simData.mpmSpacing;
    nGrad[2] = 0.25 * nzI * (1+nxI*natCod.x) * (1+nyI*natCod.y) / simData.mpmSpacing;

    //val = 0.125 * (1+natCod.x*nxI) * (1+natCod.y*nyI) * (1+natCod.z*nzI);
}

__device__ void InterpolateParticleStrain(cfloat3& xp, uint i, bufList buf,float* vgrad){
    int mx,my,mz;
    mx = (xp.x - simData.minVec.x)/simData.mpmSpacing;
    my = (xp.y - simData.minVec.y)/simData.mpmSpacing;
    mz = (xp.z - simData.minVec.z)/simData.mpmSpacing;
    int xl = simData.mpmXl;
    int yl = simData.mpmYl;
    int zl = simData.mpmZl;

    int xid = my * xl * zl + mz * xl + mx; //(-1,-1,-1 grid point)
    cfloat3 xI = buf.mpmPos[xid];
    cfloat3 center = cfloat3(xI.x+simData.mpmSpacing*0.5, xI.y+simData.mpmSpacing*0.5, xI.z+simData.mpmSpacing*0.5);
    cfloat3 natCod = (xp-center)*2/simData.mpmSpacing;

    float dxN[3];
    int tmpid;
    float vI[3];
    float tmpvgrad[9];

    for (int nx = -1; nx<=1; nx+=2){
        for (int ny=-1; ny<=1; ny+=2){
            for (int nz=-1; nz<=1; nz+=2){
                tmpid = xid;
                if (nx>0)
                    tmpid += 1;
                if (ny>0)
                    tmpid += xl*zl;
                if (nz>0)
                    tmpid += xl;

                vI[0] = buf.mpmVel[tmpid*2].x;
                vI[1] = buf.mpmVel[tmpid*2].y;
                vI[2] = buf.mpmVel[tmpid*2].z;

                dxN[0] = 0.25 * nx * (1+ny*natCod.y) * (1+nz*natCod.z) / simData.mpmSpacing;
                dxN[1] = 0.25 * ny * (1+nx*natCod.x) * (1+nz*natCod.z) / simData.mpmSpacing;
                dxN[2] = 0.25 * nz * (1+nx*natCod.x) * (1+ny*natCod.y) / simData.mpmSpacing;

                TensorProd(vI, dxN, tmpvgrad);
                for (int mi=0; mi<9; mi++){
                    vgrad[mi] += tmpvgrad[mi];
                }
            }
        }
    }
}

__device__ void InterpolateParticlePosVel(cfloat3& xp, uint i, bufList buf, cfloat3& vel, cfloat3& accel){
    int mx,my,mz;
    mx = (xp.x - simData.minVec.x)/simData.mpmSpacing;
    my = (xp.y - simData.minVec.y)/simData.mpmSpacing;
    mz = (xp.z - simData.minVec.z)/simData.mpmSpacing;
    int xl = simData.mpmXl;
    int yl = simData.mpmYl;
    int zl = simData.mpmZl;

    int xid = my * xl * zl + mz * xl + mx; //(-1,-1,-1 grid point)
    cfloat3 xI = buf.mpmPos[xid];
    cfloat3 center = cfloat3(xI.x+simData.mpmSpacing*0.5, xI.y+simData.mpmSpacing*0.5, xI.z+simData.mpmSpacing*0.5);
    cfloat3 natCod = (xp-center)*2/simData.mpmSpacing;

    float weight;
    int tmpid;
    cfloat3 velI;
    cfloat3 forceI;
    
    for (int nx = -1; nx<=1; nx+=2){
        for (int ny=-1; ny<=1; ny+=2){
            for (int nz=-1; nz<=1; nz+=2){
                tmpid = xid;
                if (nx>0)
                    tmpid += 1;
                if (ny>0)
                    tmpid += xl*zl;
                if (nz>0)
                    tmpid += xl;

                if (buf.MFtype[i]==3){ //fluid
                    velI = buf.mpmVel[tmpid*2+1];
                    forceI = buf.mpmForce[tmpid*2+1];
                }
                else{ //solid
                    velI = buf.mpmVel[tmpid*2];
                    forceI = buf.mpmForce[tmpid*2];
                }

                
                weight = 0.125 * (1+nx*natCod.x) * (1+ny*natCod.y) * (1+nz*natCod.z);
                //printf("%f %f %f    %f %f %f    %f\n",buf.mpmPos[tmpid].x,buf.mpmPos[tmpid].y,buf.mpmPos[tmpid].z, buf.mpmVel[tmpid].x,buf.mpmVel[tmpid].y,buf.mpmVel[tmpid].z,weight);
                vel += weight*velI;
                accel += weight*forceI;
            }
        }
    }
}





__device__ void contributeGridMassVel_ShapeFunction(float* mass, cfloat3* vel, float* alpha, int i, cfloat3 p, int cell, bufList buf){
    float dist[3];
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2 / d2;
	float s = 1;

	if ( buf.mgridcnt[cell] == 0 ) 
        return;
	
	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];
	float weight;

	for ( int cndx = cfirst; cndx < clast; cndx++ ) {
		uint j = buf.mgrid[cndx]; //particle index
		
        ShapeFunctionValue(p, buf.mpos[j], weight);

        
        if (buf.MFtype[j]==3){//fluid
            vel[1] += weight * buf.mf_restmass[j] * buf.mvel[j];
            mass[1] += weight* buf.mf_restmass[j];
        }
        else{//solid
            vel[0] += weight * buf.mf_restmass[j] * buf.mvel[j];
            mass[0] += weight* buf.mf_restmass[j];
        }

        for (int ki=0; ki<MAX_FLUIDNUM; ki++){
            alpha[ki] += weight* buf.mf_alpha[j*MAX_FLUIDNUM+ki];
        }
	}
}

__global__ void GetGridMassVel(bufList buf, int mpmSize){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= mpmSize ) return;

    // Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
    int gc = buf.mpmGid[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	cfloat3 pos = buf.mpmPos[ i ];
    float mass[2] = {0, 0};
    cfloat3 vel[2];
    vel[0] = make_cfloat3(0,0,0); //solid
    vel[1] = make_cfloat3(0,0,0); //fluid

    float alpha[MAX_FLUIDNUM];
    for (int ki=0; ki<MAX_FLUIDNUM; ki++)
        alpha[ki] = 0;

	for (int c=0; c < simData.gridAdjCnt; c++) {
	    contributeGridMassVel_ShapeFunction ( mass, vel, alpha, i, pos, gc + simData.gridAdj[c], buf );
	}
    buf.mpmMass[i*2] = mass[0];
    buf.mpmMass[i*2+1] = mass[1];

    if (mass[0] > 0)
        buf.mpmVel[i*2]   = vel[0]/mass[0];
	else
		buf.mpmVel[i*2] = make_cfloat3(0,0,0);
    if (mass[1] > 0)
        buf.mpmVel[i*2+1] = vel[1]/mass[1];
	else
		buf.mpmVel[i*2+1] = make_cfloat3(0,0,0);
    
    float alphasum=0;
    for (int ki=0; ki<MAX_FLUIDNUM; ki++){
        alphasum+=alpha[ki];
    }
    if (alphasum==0)
        alphasum = 1.0f;
    for (int ki=0; ki<MAX_FLUIDNUM; ki++){
        alpha[ki] /= alphasum;
        buf.mpmAlpha[i*MAX_FLUIDNUM+ki] = alpha[ki];
    }
}

/*
__global__ void CalcMpmParticleTensor(bufList buf, int pnum){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;
    if (buf.MFtype[i]!=1) return;

    // Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
    int gc = buf.mgcell[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

    //Velocity Gradient
    float velGrad[9], velGradT[9];
    for (int mi=0; mi<9; mi++)
        velGrad[mi] = 0;
    
    InterpolateParticleStrain(buf.mpos[i], i, buf, velGrad);
	
    transpose3(velGrad, velGradT);
    
    //Calculate Tensor
    float tensorS[9], tensorFull[9];
    float epsilon[9], omega[9];
    float dTensorS[9];
    for (int mi=0; mi<9; mi++)
    {
        tensorS[mi] = buf.MFtensor[i*9+mi];
        //tensorFull[mi] = buf.MFtensor[i*9+mi];
        epsilon[mi] = (velGrad[mi] + velGradT[mi])*0.5;
        omega[mi] = (velGrad[mi] - velGradT[mi])*0.5;
        buf.MFvelgrad[i*9+mi] = velGrad[mi]; //Store the VelGrad Value
        dTensorS[mi] = 0;
    }
    
    //buf.mclr[i] = COLORA(abs(velGrad[0]),abs(velGrad[4]),abs(velGrad[8]),1);

    float coG = simData.coG;
    float coK = simData.coK;
    //update tensor
    for (int ii=0; ii<3; ii++){
        for (int ij=0; ij<3; ij++){
            for (int ik=0; ik<3; ik++){
                dTensorS[ii*3+ij] += tensorS[ii*3+ik]*omega[ij*3+ik] + tensorS[ij*3+ik]*omega[ii*3+ik];
            }
        }
    }
    
    float epsum = (epsilon[0] + epsilon[4] + epsilon[8])/3.0;
    
    epsilon[0] -= epsum;
    epsilon[4] -= epsum;
    epsilon[8] -= epsum;
    for (int mi=0; mi<9; mi++){
        dTensorS[mi] += 2 * coG * epsilon[mi];
        buf.MFtensor[i*9+mi] += dTensorS[mi] * simData.mf_dt;
    }
    epsum = epsum * coK * simData.mf_dt;
    buf.MFtensor[i*9]   += epsum;
    buf.MFtensor[i*9+4] += epsum;
    buf.MFtensor[i*9+8] += epsum;

	//float checkcolor = abs(buf.MFtensor[i*9+0])+abs(buf.MFtensor[i*9+1])+abs(buf.MFtensor[i*9+2]);
	//buf.mclr[i] = COLORA( abs(buf.MFtensor[i*9+0])/checkcolor,abs(buf.MFtensor[i*9+1])/checkcolor,abs(buf.MFtensor[i*9+2])/checkcolor,1);


}

__device__ void mul_mat_vec(float* a, float* b, float* c){
    //a:mat; b:vec; c:vec
    float tmp[3];
    tmp[0] = a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
    tmp[1] = a[3]*b[0] + a[4]*b[1] + a[5]*b[2];
    tmp[2] = a[6]*b[0] + a[7]*b[1] + a[8]*b[2];
    c[0] = tmp[0];
    c[1] = tmp[1];
    c[2] = tmp[2];
}



__device__ void GridCollisionDetection(cfloat3& pos, cfloat3& vel){
    // Get particle vars
	register cfloat3 accel, norm;
	register float diff, adj, speed;
    float time = 0;
    float ss = simData.psimscale;

	accel = make_cfloat3(0,0,0);

	// Soft Boundaries
	// Y-axis
	diff = simData.pradius - (pos.y - (simData.pboundmin.y + (pos.x-simData.pboundmin.x)*simData.pground_slope )) * ss;

	if ( diff > EPSILON ) {
		norm = make_cfloat3( -simData.pground_slope, 1.0 - simData.pground_slope, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;
	}

	diff = simData.pradius - ( simData.pboundmax.y - pos.y )*ss;

	if ( diff > EPSILON ) {
		norm = make_cfloat3(0, -1, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;
	}

	
#ifdef _xzsoftmargin
	// X-axis
	diff = simData.pradius - (pos.x - (simData.pboundmin.x + (sin(time*simData.pforce_freq)+1)*0.5 * simData.pforce_min))*ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3( 1, 0, 0);
		adj = (simData.pforce_min+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}
	diff = simData.pradius - ( (simData.pboundmax.x - (sin(time*simData.pforce_freq)+1)*0.5*simData.pforce_max) - pos.x)*ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3(-1, 0, 0);
		adj = (simData.pforce_max+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}

	// Z-axis
	diff = simData.pradius - (pos.z - simData.pboundmin.z ) * ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3( 0, 0, 1 );
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}
	diff = simData.pradius - ( simData.pboundmax.z - pos.z )*ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3( 0, 0, -1 );
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}
#endif
	
	// Leap-frog Integration
	vel = accel*simData.mf_dt + vel;				// v(t+1/2) = v(t-1/2) + a(t) dt
}

__device__ void ContributeMpmGridForce_ShapeFunction(cfloat3* res, uint i, uint gcell, bufList buf){
    cfloat3 p = buf.mpmPos[i];

    float dist[3], absdist[3];
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2 / d2;
	float s = 1;

	int cfirst = buf.mgridoff[ gcell ];
	int clast = cfirst + buf.mgridcnt[ gcell ];
	float weight;
    float DN[3];
    float tmpMat[9];
    float tensorF[3]={0,0,0};
    cfloat3 pressureForce = make_cfloat3(0,0,0);
    cfloat3 viscForce = make_cfloat3(0,0,0);
    cfloat3 bodyForce0 = make_cfloat3(0,0,0);
    cfloat3 bodyForce1 = make_cfloat3(0,0,0);
    cfloat3 deformForce = make_cfloat3(0,0,0);
    float vj;

    bool isLiquid = (buf.mpmAlpha[i*MAX_FLUIDNUM+1]<0.5);

	for ( int cndx = cfirst; cndx < clast; cndx++ ) {
		uint j = buf.mgrid[cndx];
        
        vj = buf.mf_restmass[j] / buf.mf_restdensity[j];

        ShapeFunctionValue(p, buf.mpos[j], weight);
        ShapeFunctionGradValue(p, buf.mpos[j], DN);

        //Add pressure force
		pressureForce.x += buf.mpress[j] * DN[0] * vj;
		pressureForce.y += buf.mpress[j] * DN[1] * vj;
		pressureForce.z += buf.mpress[j] * DN[2] * vj;
		

        //Add viscosity force
        viscForce.x += (-1) * buf.mf_visc[j] * DN[0] * buf.MFvelgrad[j*9] * vj;
        viscForce.y += (-1) * buf.mf_visc[j] * DN[1] * buf.MFvelgrad[j*9+4] * vj;
        viscForce.z += (-1) * buf.mf_visc[j] * DN[2] * buf.MFvelgrad[j*9+8] * vj;

        //Tensor Force
        if (buf.MFtype[j]==1){
            mul_mat_vec( &buf.MFtensor[j*9], DN, tensorF);
            deformForce.x -=  tensorF[0] * vj;
            deformForce.y -=  tensorF[1] * vj;
            deformForce.z -=  tensorF[2] * vj;
        }
        
        //body Force
        if (buf.MFtype[j]==1)//solid
            bodyForce0 += buf.mf_restmass[j] * (simData.pgravity) * weight;
        else
            bodyForce1 += buf.mf_restmass[j] * simData.pgravity * weight;
            
    }

    res[1] += pressureForce;
    res[1] += viscForce;

    res[0] += deformForce;
    
    res[0] += bodyForce0;
    res[1] += bodyForce1;
}

__device__ void ContributeMpmGridMassGrad(uint i, uint cell, bufList buf, cfloat3& massGrad){
	cfloat3 p = buf.mpmPos[i];

	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];
    float DN[3];
    
	for ( int cndx = cfirst; cndx < clast; cndx++ ) {
		uint j = buf.mgrid[cndx];
		if(buf.MFtype[j]!=1)
			continue;

        ShapeFunctionGradValue(p, buf.mpos[j], DN);
        //Add pressure force
        massGrad.x += (-1)*buf.mf_restmass[j] * DN[0];
        massGrad.y += (-1)*buf.mf_restmass[j] * DN[1];
        massGrad.z += (-1)*buf.mf_restmass[j] * DN[2];
	}
}

//Make the velocity converge
__device__ void FluidSolidContact(uint i, bufList buf, cfloat3& force){
	//Judge If Two Phases Meet
	if(buf.mpmMass[i*2] == 0 || buf.mpmMass[i*2+1]==0)
		return;

	//Compute Mass Gradient = Normal Direction
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
    int gc = buf.mpmGid[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	cfloat3 massGrad = make_cfloat3(0,0,0);
	for (int c=0; c < simData.gridAdjCnt; c++) {
	    ContributeMpmGridMassGrad ( i, gc + simData.gridAdj[c], buf, massGrad );
	}

	//Judge If Collision Happens
	cfloat3 dv = buf.mpmVel[2*i+1] - buf.mpmVel[2*i];
	float dvMG = dot(dv, -massGrad);
	if(dvMG >= 0)
		return;

	//Compute Force For Liquid ( Minus Solid Force )
	force = -dv * (buf.mpmMass[2*i]*buf.mpmMass[2*i+1])/(buf.mpmMass[2*i]+buf.mpmMass[2*i+1])/simData.mf_dt;

}


__global__ void CalcMpmGridForce(bufList buf, int mpmSize){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= mpmSize ) return;

    float nodemass = buf.mpmMass[i*2] + buf.mpmMass[i*2+1];
    if (nodemass == 0)
        return;

    // Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
    int gc = buf.mpmGid[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

    //Velocity Gradient
    cfloat3 gForce[2];
    gForce[0] = make_cfloat3(0,0,0);
    gForce[1] = make_cfloat3(0,0,0);
    
    bool isLiquid = (buf.mpmAlpha[i*MAX_FLUIDNUM+1]<0.5);
	
    for (int c=0; c < simData.gridAdjCnt; c++) {
	    ContributeMpmGridForce_ShapeFunction ( gForce, i, gc + simData.gridAdj[c], buf );
	}

    buf.mpmForce[i*2] = gForce[0] / buf.mpmMass[i*2];
    if (buf.mpmMass[i*2]==0)
        buf.mpmForce[i*2] = make_cfloat3(0,0,0);

    buf.mpmForce[i*2+1] = gForce[1] / buf.mpmMass[i*2+1];
    if (buf.mpmMass[i*2+1]==0)
        buf.mpmForce[i*2+1] = make_cfloat3(0,0,0);

	cfloat3 contactLiquidForce = make_cfloat3(0,0,0);
	FluidSolidContact(i, buf, contactLiquidForce); //For Now No Friction
	buf.mpmForce[i*2] -= contactLiquidForce/buf.mpmMass[2*i] * paramCarrier.collisionStiff;
	buf.mpmForce[i*2+1] += contactLiquidForce/buf.mpmMass[2*i+1] * paramCarrier.collisionStiff;


    buf.mpmVel[i*2] = buf.mpmVel[i*2] + buf.mpmForce[i*2] * simData.mf_dt;
    buf.mpmVel[i*2+1] = buf.mpmVel[i*2+1] + buf.mpmForce[i*2+1] * simData.mf_dt;

    //GridCollisionDetection(buf.mpmPos[i], buf.mpmVel[i]);
    return;
    
}

__device__ void CalcMpmParticleCollision(bufList buf, uint i, cfloat3& accel, cfloat3 vel){
    
	register cfloat3 pos = buf.mpos[i];
	accel = make_cfloat3(0,0,0);
    
    float dx;
    float stiff = simData.pbstiff;
    float damp = simData.pdamp;
	float margin = 0;

    if (buf.misbound[i] != 1)
	{
        if (buf.mpos[i].y < simData.pboundmin.y + margin && vel.y<0){
            dx = simData.pboundmin.y + margin - buf.mpos[i].y;
            accel += (stiff * dx  - damp*vel.y)* make_cfloat3(0,1,0);
        }

        if (buf.mpos[i].y > simData.pboundmax.y - margin && vel.y>0){
            dx = buf.mpos[i].y - simData.pboundmax.y + margin ;
            accel += (stiff * dx + damp * vel.y) * make_cfloat3(0,-1,0);
        }

        if (buf.mpos[i].x < simData.pboundmin.x + margin && vel.x<0){
            dx = simData.pboundmin.x + margin - buf.mpos[i].x;
            accel += stiff * dx * make_cfloat3(1,0,0);
        }
        
        if (buf.mpos[i].x > simData.pboundmax.x - margin && vel.x>0){
            dx = buf.mpos[i].x - simData.pboundmax.x + margin;
            accel += stiff * dx * make_cfloat3(-1,0,0);
        }

        if (buf.mpos[i].z < simData.pboundmin.z + margin && vel.z<0){
            dx = simData.pboundmin.z + margin - buf.mpos[i].z;
            accel += stiff * dx * make_cfloat3(0,0,1);
        }
        
        if (buf.mpos[i].z > simData.pboundmax.z - margin && vel.z>0){
            dx = buf.mpos[i].z - simData.pboundmax.z + margin;
            accel += stiff * dx * make_cfloat3(0,0,-1);
        }
	}
}


__device__ void CollisionDetection ( cfloat3& pos, cfloat3& vel, cfloat3& accel0 )
{		
	// Get particle vars
	register cfloat3 accel, norm;
	register float diff, adj, speed;
	register float newdens,newvisc;
    float time = 0;
    float ss = simData.psimscale;

	accel = make_cfloat3(0,0,0);

	// Soft Boundaries
	// Y-axis
	diff = simData.pradius - (pos.y - (simData.pboundmin.y + (pos.x-simData.pboundmin.x)*simData.pground_slope )) * ss;

	if ( diff > EPSILON ) {
		norm = make_cfloat3( -simData.pground_slope, 1.0 - simData.pground_slope, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;
	}

	diff = simData.pradius - ( simData.pboundmax.y - pos.y )*ss;

	if ( diff > EPSILON ) {
		norm = make_cfloat3(0, -1, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;
	}

	
#ifdef _xzsoftmargin
	// X-axis
	diff = simData.pradius - (pos.x - (simData.pboundmin.x + (sin(time*simData.pforce_freq)+1)*0.5 * simData.pforce_min))*ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3( 1, 0, 0);
		adj = (simData.pforce_min+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}
	diff = simData.pradius - ( (simData.pboundmax.x - (sin(time*simData.pforce_freq)+1)*0.5*simData.pforce_max) - pos.x)*ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3(-1, 0, 0);
		adj = (simData.pforce_max+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}

	// Z-axis
	diff = simData.pradius - (pos.z - simData.pboundmin.z ) * ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3( 0, 0, 1 );
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}
	diff = simData.pradius - ( simData.pboundmax.z - pos.z )*ss;
//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if ( diff > EPSILON ) {
		norm = make_cfloat3( 0, 0, -1 );
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, vel );
		norm *= adj; accel += norm;//*scale_dens;
	}
#endif
	
	//accel += accel0;
	accel0 += accel;



	// Leap-frog Integration

	vel = accel*simData.mf_dt + vel;				// v(t+1/2) = v(t-1/2) + a(t) dt
}



__global__ void UpdateMpmParticlePos(bufList buf, int pnum){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

    // Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
    int gc = buf.mgcell[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

    cfloat3 vel = make_cfloat3(0,0,0);
    cfloat3 accel = make_cfloat3(0,0,0);

    InterpolateParticlePosVel(buf.mpos[i], i, buf, vel, accel);

    CollisionDetection(buf.mpos[i], vel, accel);
    
	cfloat3 vpic,vflip;
    vpic = vel; //grid interpolation
    vflip = buf.mvel[i] + accel * simData.mf_dt;
	float a;
    if(buf.MFtype[i]==3) 
        a = paramCarrier.LFlipControl;
    else
        a = paramCarrier.SFlipControl;
    
    buf.mvel[i] = a*vpic + (1-a)*vflip;
    buf.mveleval[i] = buf.mvel[i];
    buf.mpos[i] += buf.mvel[i] * simData.mf_dt / simData.psimscale;
	buf.mforce[i] = simData.pgravity - accel; //for drift velocity

}

*/