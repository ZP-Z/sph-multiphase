
#include "fluid_system_kern.cuh"
#include "Mobility.cuh"
#include <stdio.h>
//#include "stokesian\NeighborGrid.cuh"

__device__ float* signature;
extern __device__ ParamCarrier paramCarrier;
extern bufList	fbuf;

__inline__ __device__ __host__ int tsId(int i, int j, int k) {
	return i*9 + j*3 + k;
}

__inline__ __device__ __host__ int kroneckerDelta(int i, int j) {
	return i==j?1:0;
}

void setSignature() {
	float* ptr;
	cudaMalloc(&ptr, 27*sizeof(float));
	cudaMemcpyToSymbol(signature, &ptr, sizeof(ptr));
	float val[27];
	for(int i=0; i<27; i++)
		val[i] = 0;
	val [tsId(0,1,2)] = 1;
	val [tsId(2,0,1)] = 1;
	val [tsId(1,2,0)] = 1;
	val [tsId(0,2,1)] = -1;
	val [tsId(2,1,0)] = -1;
	val [tsId(1,0,2)] = -1;
	cudaMemcpy(ptr, val, sizeof(float)*27, cudaMemcpyHostToDevice);
	return;
}

__device__ void PairwiseMobMatrix(cfloat3 drv, float A11[],float A1N[], float C11[], float C1N[], float B11[], float B1N[]) {
	float dr = sqrt(dot(drv, drv));
	float eu[3];
	eu[0] = drv.x/dr;
	eu[1] = drv.y/dr;
	eu[2] = drv.z/dr;

	//get XYABC mobility functions
	float X11A = 1;
	float dr3 = dr*dr*dr;
	float X1NA = 1.5/dr - 1/dr3;
	float Y11A = 1;
	float Y1NA = 0.75/dr + 0.5/dr3;
	float Y11B = 0;
	float Y1NB = -3/(4*dr*dr);
	float X11C = 0.75;
	float X1NC = 0.75/dr3;
	float Y11C = 0.75;
	float Y1NC = -0.375/dr3;

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			A11[i*3+j] = X11A * eu[i]*eu[j] + Y11A * (kroneckerDelta(i, j)-eu[i]*eu[j]);
			A1N[i*3+j] = X1NA * eu[i]*eu[j] + Y1NA * (kroneckerDelta(i, j)-eu[i]*eu[j]);
			C11[i*3+j] = X11C * eu[i]*eu[j] + Y11C * (kroneckerDelta(i, j)-eu[i]*eu[j]);
			C1N[i*3+j] = X1NC * eu[i]*eu[j] + Y1NC * (kroneckerDelta(i, j)-eu[i]*eu[j]);
		}
	}

	for (int i=0; i<3; i++) {
		A11[i*3+i] -= 1;
		C11[i*3+i] -= 0.75;
	}

	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			float dot = 0;
			for (int k=0; k<3; k++) {
				dot = dot + signature[tsId(j, i, k)]*eu[k];
			}
			B11[i*3+j] = Y11B * dot;
			B1N[i*3+j] = -Y1NB * dot;
		}
	}


}


__global__ void MVmultiply(float* dst, float* mat, float* f,int n) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	int Adim = 6*n;
	if (i>=Adim)
		return;
	
	float sum = 0;
	for (int j=0; j<Adim; j++) {
		sum += mat[i*Adim+j] * f[j];
	}
	dst[i] = sum;
	//printf("%d %f %f\n",i, f[i], dst[i]);
	return;
}


//calculate velocity
void getMobU(float* mat,float* f, float* u, int nsize) {
	//simply matrix multiply
	int threadnum = 256;
	int blocknum = (6*nsize-1)/threadnum + 1;

	MVmultiply<<<blocknum, threadnum>>>(u, mat, f, nsize);
	cudaDeviceSynchronize();
}



//-------------------------------------------------
//-------------------------------------------------
//-------------                   -----------------
//-------------   cutoff version  -----------------
//-------------                   -----------------
//-------------------------------------------------
//-------------------------------------------------

//math functions
__inline__ __device__ __host__ void mulmv3(float* dst, float* mat, float* v) {
	float tmp[3] ={0,0,0};
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			tmp[i] += mat[i*3+j]*v[j];
		}
	}
	dst[0] = tmp[0];
	dst[1] = tmp[1];
	dst[2] = tmp[2];
}

__device__ __host__ void mulmat(float* dst, float* a,float* b){
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			dst[i*3+j]=0;
			for(int k=0;k<3;k++){
				dst[i*3+j]+=a[i*3+k]*b[k*3+j];
			}
		}
	}
}

//__global__ void Kern_getmobu_cutoff(float* pos, int nsize, float* f, float* u, sortingGrid grid) {
//	int i = blockIdx.x * blockDim.x + threadIdx.x;
//	if(i>=nsize)
//		return;
//
//	//locate i's cell
//	int id = grid.idsort[i];
//	int cell = grid.cellId[i];
//	if(cell==GRIDUNDEF)
//		return;
//
//	float3 ipos = grid.pos[id];
//	float3 dx;
//	int ncell;
//
//	int kid;
//	float3 drv;
//	int n1,n2,n3,n4;
//	float mat[9];
//	memset(mat,0,9*sizeof(float));
//	float3 utmp = make_float3(0,0,0);
//	float unity[9] = {1,0,0, 0,1,0, 0,0,1};
//	//loop through neighboring cells
//	int count = 0;
//	u[3*id] = 0;
//	u[3*id+1] = 0;
//	u[3*id+2] = 0;
//	
//	for (int j=0; j<27; j++) {
//		ncell = cell + grid.searchList[j];
//
//		count += grid.gridend[ncell] - grid.gridstart[ncell];
//		for (int k=grid.gridstart[ncell]; k<grid.gridend[ncell]; k++) {
//			kid = grid.idsort[k];
//			if(kid != id){
//				drv = grid.pos[kid] - ipos;
//				float A11[9], A1N[9], C11[9], C1N[9], B11[9], B1N[9];
//
//				PairwiseMobMatrix(drv, A11, A1N, C11, C1N, B11, B1N);
//				mulmv3((float*)&utmp, A1N, f+3*kid);
//				u[3*id] += utmp.x;
//				u[3*id+1] += utmp.y;
//				u[3*id+2] += utmp.z;
//
//				/*if (id==0) {
//					printf("%d %d %f %f %f\n",id, kid, grid.pos[kid].x, grid.pos[kid].y, grid.pos[kid].z);
//				}*/
//			}
//			else {
//				mulmv3((float*)&utmp, unity, f+3*kid);
//				u[3*id] += utmp.x;
//				u[3*id+1] += utmp.y;
//				u[3*id+2] += utmp.z;
//			}
//		}
//		
//	}
//	//printf("%d %d\n", id, count);
//}



__global__ void Kern_getmobu_walkthrough( stokesianBufList buflist, int nsize) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i>=nsize)
		return;
	
	int id = i;
	cfloat3 ipos = buflist.dispBuffer[id].pos;
	int kid;
	cfloat3 drv;

	cfloat3 utmp(0, 0, 0);
	float unity[9] ={1,0,0, 0,1,0, 0,0,1};
	float force[3];
	

	//buflist.cuUnew[id].Set(0,0,0);
	buflist.cuOmega[id].Set(0,0,0);
	
	
	for (int k=0; k<nsize; k++) {
		kid = k;
		if (kid != id) {
			drv = buflist.dispBuffer[kid].pos - ipos;
			
			
			float A11[9], A1N[9], C11[9], C1N[9], B11[9], B1N[9];

			PairwiseMobMatrix(drv, A11, A1N, C11, C1N, B11, B1N);

			//get velocity
			force[0] = buflist.cuForce[kid].x;
			force[1] = buflist.cuForce[kid].y;
			force[2] = buflist.cuForce[kid].z;

			mulmv3((float*)&utmp, A1N, force);
			buflist.cuUnew[id] += utmp;

			//get angle velocity
			mulmv3((float*)&utmp, B1N, force);
			buflist.cuOmega[id] += utmp;
		}
		else {
			mulmv3((float*)&utmp, unity, (float*) &buflist.cuForce[kid]);
			buflist.cuUnew[id] += utmp;
		}
	}
	
	buflist.dispBuffer[id].pos += buflist.cuUnew[id] * buflist.dt;
}

__global__ void Kern_getForce(stokesianBufList buflist, int nsize){
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i>=nsize)
		return;

	int id = i;
	cfloat3 ipos = buflist.dispBuffer[id].pos;
	int kid;
	cfloat3 drv;

	
	int count = 0;
	for (int k=0; k<nsize; k++) {
		kid = k;
		if (kid != id) {
			drv = buflist.dispBuffer[kid].pos - ipos;

			//near-field interactions between particles
			float dr = sqrt(dot(drv, drv));
			
			//splitting
			
			float actdist = 1.5;
			if (dr < actdist) {
				float fac = 5*(exp(actdist/dr)-2.7128*dr/actdist);
				cfloat3 intforce = drv * (-1) * fac;
				buflist.cuForce[i] += intforce;
			}

			//surface tension
			//float actdist2 = 4;
			//if (dr<actdist2) {
			//	float fac = 10*cos(1.5*3.14159*dr/actdist2);
			//	if(fac>0)
			//		fac *= 30;
			//	//printf("%f %f %f %f %f %f\n",drv.x, drv.y, drv.z, dr, fac);
			//	float3 intforce = -fac * drv/dr/dr;
			//	f[i*3+0] += intforce.x;
			//	f[i*3+1] += intforce.y;
			//	f[i*3+2] += intforce.z;
			//}
		}
	}
}

//void getMobU_cutoff(float* pos, int nsize, float* f, float* u, sortingGrid grid) {
//	int threadnum = 256;
//	int blocknum = (nsize-1)/256+1;
//
//	//Kern_getmobu_cutoff<<<blocknum,threadnum>>>(pos, nsize, f, u, grid);
//	Kern_getmobu_walkthrough<<<blocknum, threadnum>>>(pos, nsize, f, u);
//	cudaDeviceSynchronize();
//	
//	return;
//}

__device__ __inline__ float getWeight(cfloat3 pos, cint3 idvec, float cellsize) {
	cfloat3 dx = pos - (cfloat3)idvec * cellsize;
	return (1 - abs(dx.x)/cellsize) * (1 - abs(dx.y)/cellsize) * (1 - abs(dx.z)/cellsize);
}

__device__ __inline__ int getId(cint3 idvec) {
	return idvec.y*paramCarrier.mpmRes.x*paramCarrier.mpmRes.z + idvec.z*paramCarrier.mpmRes.x + idvec.x;
}

__global__ void Kern_interpolateVelocity(stokesianBufList buflist, int pnum, bufList fbuf) {
	int i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i>=pnum)
		return;

	cfloat3 ipos = buflist.dispBuffer[i].pos;
	cint3 mpmres = paramCarrier.mpmRes;
	float cellsize = paramCarrier.mpmcellsize;
	
	ipos = ipos - paramCarrier.mpmXmin;
	cfloat3 idvecf = ipos / paramCarrier.mpmcellsize;
	cint3 idvec(idvecf.x, idvecf.y, idvecf.z);
	cfloat3 dx;
	float weight;
	cfloat3 velsum(0,0,0);

	
	cint3 tmp3;
	tmp3 = idvec + cint3(0, 0, 0);
	velsum  += fbuf.mpmVel[ getId(tmp3) ] * getWeight(ipos, tmp3, cellsize);
	
	tmp3 = idvec + cint3(1, 0, 0);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	tmp3 = idvec + cint3(0, 1, 0);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	tmp3 = idvec + cint3(0, 0, 1);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	tmp3 = idvec + cint3(1, 1, 0);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	tmp3 = idvec + cint3(1, 0, 1);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	tmp3 = idvec + cint3(0, 1, 1);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	tmp3 = idvec + cint3(1, 1, 1);
	velsum  += fbuf.mpmVel[getId(tmp3)] * getWeight(ipos, tmp3, cellsize);

	//printf("%f %f %f\n",velsum.x, velsum.y,velsum.z);
	buflist.cuUnew[i] = velsum / paramCarrier.simscale;
}


void getMobU_walkthrough(stokesianBufList buflist, int pnum) {
	int threadnum = 256;
	int blocknum = (pnum-1)/256+1;

	Kern_getForce<<<blocknum, threadnum>>>(buflist, pnum);
	cudaDeviceSynchronize();

	Kern_interpolateVelocity<<<blocknum, threadnum>>>(buflist, pnum, fbuf);
	cudaDeviceSynchronize();

	Kern_getmobu_walkthrough<<<blocknum, threadnum>>>(buflist, pnum);
	cudaDeviceSynchronize();

	return;
}