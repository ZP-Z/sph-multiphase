

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <conio.h>


#include "fluid_system_host.cuh"		
#include "fluid_system_kern.cuh"
#include "MpmSolver.cuh"
#include "sph_solid.cuh"

#include "thrust\device_vector.h"	//thrust libs
#include "thrust\sort.h" 


FluidParams		fcuda;
bufList			fbuf;

cudaError_t error;
extern ParamCarrier hostCarrier;

void cudaExit (int argc, char **argv)
{
	//CUT_EXIT(argc, argv); 
}
void cudaInit(int argc, char **argv)
{   
	//CUT_DEVICE_INIT(argc, argv);
	
	cudaDeviceProp p;
	cudaGetDeviceProperties ( &p, 0);
	
	printf ( "-- CUDA --\n" );
	printf ( "Name:       %s\n", p.name );
	printf ( "Revision:   %d.%d\n", p.major, p.minor );
	printf ( "Global Mem: %d\n", p.totalGlobalMem );
	printf ( "Shared/Blk: %d\n", p.sharedMemPerBlock );
	printf ( "Regs/Blk:   %d\n", p.regsPerBlock );
	printf ( "Warp Size:  %d\n", p.warpSize );
	printf ( "Mem Pitch:  %d\n", p.memPitch );
	printf ( "Thrds/Blk:  %d\n", p.maxThreadsPerBlock );
	printf ( "Const Mem:  %d\n", p.totalConstMem );
	printf ( "Clock Rate: %d\n", p.clockRate );	
};

int iDivUp (int totalnum, int threadnum) {
	
	if(threadnum==0)
		return 1;

	return (totalnum % threadnum != 0) ? (totalnum / threadnum + 1) : (totalnum / threadnum);
}

inline bool isPowerOfTwo(int n) { return ((n&(n-1))==0) ; }
inline int floorPow2(int n) {
	#ifdef WIN32
		return 1 << (int)logb((float)n);
	#else
		int exp;
		frexp((float)n, &exp);
		return 1 << (exp - 1);
	#endif
}

// Compute number of blocks to create
void computeNumBlocks (int numPnts, int maxThreads, int &numBlocks, int &numThreads)
{
	numThreads = min( maxThreads, numPnts );
	numBlocks = iDivUp ( numPnts, numThreads );
	if(numThreads==0)
		numThreads = 1;
}
#define CUDA_SAFE_CALL

void FluidClearCUDA ()
{
	cudaFree(fbuf.displayBuffer);
	cudaFree(fbuf.calcBuffer);

	cudaFree ( fbuf.msortbuf );	
	cudaFree(fbuf.MFidTable);
	//new sort
	cudaFree(fbuf.mgcell);
	cudaFree(fbuf.mgndx);
	cudaFree(fbuf.mgridcnt);
	cudaFree ( fbuf.midsort );
	cudaFree ( fbuf.mgridoff );

}









void FluidSetupCUDA(ParamCarrier& params){
	fcuda.pnum = params.num;
	fcuda.gridTotal = params.gridtotal;
	
	// Compute number of blocks and threads
	computeNumBlocks(fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);			// particles
	computeNumBlocks(fcuda.gridTotal, 384, fcuda.gridBlocks, fcuda.gridThreads);		// grid cell

	// Allocate particle buffers
	//fcuda.szPnts = (fcuda.numBlocks  * fcuda.numThreads);
	fcuda.szPnts = params.maxNum;

	cudaMalloc(&fbuf.displayBuffer, fcuda.szPnts * sizeof(displayPack));
	cudaMalloc(&fbuf.calcBuffer,	fcuda.szPnts*sizeof(calculationPack));
	int temp_size = (sizeof(displayPack) + sizeof(calculationPack));
	cudaMalloc(&fbuf.msortbuf, fcuda.szPnts*temp_size);

	//without index sort
	cudaMalloc(&fbuf.densityResidue, fcuda.szPnts * sizeof(float));
	cudaMalloc(&fbuf.press_l,  fcuda.szPnts*sizeof(float));
	cudaMalloc(&fbuf.press_l1, fcuda.szPnts*sizeof(float));
	cudaMalloc(&fbuf.rho_adv, fcuda.szPnts*sizeof(float));
	cudaMalloc(&fbuf.stress,  fcuda.szPnts*sizeof(cmat3));
	cudaMalloc(&fbuf.aii, fcuda.szPnts*sizeof(float));
	cudaMalloc(&fbuf.dii, fcuda.szPnts*sizeof(cfloat3));
	cudaMalloc(&fbuf.dijpj, fcuda.szPnts*sizeof(cfloat3));

	cudaMalloc(&fbuf.MFidTable,	fcuda.szPnts*sizeof(int));
	

	// Allocate grid
	fcuda.szGrid = (fcuda.gridBlocks * fcuda.gridThreads);
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgcell, fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgndx,  fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridcnt, fcuda.szGrid*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.midsort,  fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridoff, fcuda.szGrid*sizeof(int)));
	

	//MpmAllocateBuffer();
}


void GetParticleIndexCUDA()
{
	cudaMemset ( fbuf.mgridcnt, 0,			fcuda.gridTotal * sizeof(int));
	cudaMemset ( fbuf.mgridoff, 0,			fcuda.gridTotal * sizeof(int));
	cudaMemset ( fbuf.mgcell, 0,			fcuda.pnum * sizeof(uint));
	InitialSort<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr,  "CUDA ERROR: InsertParticlesCUDA: %s\n", cudaGetErrorString(error) );
	}  
	cudaThreadSynchronize ();

}
void GetGridListCUDA()
{

	thrust::device_ptr<uint> dev_keysg(fbuf.mgcell);
	thrust::device_ptr<uint> dev_valuesg(fbuf.midsort);
	thrust::sort_by_key(dev_keysg,dev_keysg+fcuda.pnum,dev_valuesg);
	//cudaThreadSynchronize ();
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf(stderr, "CUDA ERROR: Thrust sort: %s\n", cudaGetErrorString(error));
	}

	CalcFirstCnt <<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	//	cudaThreadSynchronize ();
	cudaThreadSynchronize ();


	GetCnt <<<fcuda.numBlocks,fcuda.numThreads>>> (fbuf,fcuda.pnum);
	cudaThreadSynchronize ();
	
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf(stderr, "CUDA ERROR: Sort Grid: %s\n", cudaGetErrorString(error));
	}
}

void RearrageDataCUDA()
{
	// Transfer particle data to temp buffers
	int n = fcuda.pnum;
	
	cudaMemcpy ( fbuf.msortbuf + n*BUF_DISPLAYBUF, fbuf.displayBuffer, n*sizeof(displayPack), cudaMemcpyDeviceToDevice);
	cudaMemcpy(fbuf.msortbuf + n*BUF_CALCBUF, fbuf.calcBuffer, n*sizeof(calculationPack), cudaMemcpyDeviceToDevice);
	//cudaMemcpy(fbuf.msortbuf + n*BUF_INTMBUF, fbuf.intmBuffer, n*sizeof(IntermediatePack), cudaMemcpyDeviceToDevice);

	// Counting Sort - pass one, determine grid counts
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR:CopyToSortBufferCUDA: %s\n", cudaGetErrorString(error) );
	} 

	RearrangeData <<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	cudaThreadSynchronize ();

	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR:Sorting Failed: %s\n", cudaGetErrorString(error) );
	} 

}






void MpmAllocateBufferCUDA(ParamCarrier& param){

    //calculate size
    fcuda.mpmSize = param.mpmNodeNum;

    //int splitnum = 2;

    computeNumBlocks ( fcuda.mpmSize, 384, fcuda.mpmBlocks, fcuda.mpmThreads);

    cudaMalloc(&fbuf.mpmMass,   fcuda.mpmSize * sizeof(float));
	cudaMalloc(&fbuf.mpmPos, fcuda.mpmSize * sizeof(cfloat3));
	cudaMalloc(&fbuf.mpmVel, fcuda.mpmSize * sizeof(cfloat3));
    //cudaMalloc(&fbuf.mpmAlpha, fcuda.mpmSize * sizeof(float) * MAX_FLUIDNUM);
    //cudaMalloc(&fbuf.mpmForce, fcuda.mpmSize * sizeof(cfloat3) * splitnum);
	cudaMalloc(&fbuf.mpmGid,		fcuda.mpmSize * sizeof(uint));
}


void IndexMPMSortCUDA(){
    initMpm <<< fcuda.mpmBlocks, fcuda.mpmThreads >>> (fbuf, fcuda.mpmSize);
    cudaThreadSynchronize();
}


void MpmColorTestCUDA() {
	MpmColorTest <<< fcuda.mpmBlocks, fcuda.mpmThreads>>>(fbuf, fcuda.mpmSize);
	cudaThreadSynchronize();
}

void MpmGetMomentumCUDA() {
	MpmGetMomentum <<< fcuda.mpmBlocks, fcuda.mpmThreads>>>(fbuf, fcuda.mpmSize);
	cudaThreadSynchronize();
}


void initSPH(float* restdensity,int* mftype){
	
    //initDensity<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	//cudaThreadSynchronize();

	/*CUDA_SAFE_CALL( cudaMemcpy( restdensity, fbuf.mf_restdensity, fcuda.pnum*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL( cudaMemcpy( mftype,      fbuf.MFtype,         fcuda.pnum*sizeof(int),   cudaMemcpyDeviceToHost));

	double sum=0;
	int cnt=0;
	for(int i=0; i<fcuda.pnum; i++){
		if( mftype[i]==1){
			sum += restdensity[i];
			cnt++;
		}
	}
	if (cnt > 0)
		sum /= cnt;
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.mf_restdensity,restdensity, fcuda.pnum*sizeof(float), cudaMemcpyHostToDevice));
	printf("average density %f\n",sum);*/
    
    //MpmSortGridCuda();

    //initMpm<<<fcuda.mpmBlocks, fcuda.mpmThreads>>>(fbuf, fcuda.mpmSize);
    //cudaThreadSynchronize();	
}


void MfComputePressureCUDA ()
{
	//ComputeBoundaryVolume <<< fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	//cudaDeviceSynchronize();

	ComputeDensityPressure <<< fcuda.numBlocks, fcuda.numThreads >>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
	
	error = cudaGetLastError();
	if(error != cudaSuccess){
		printf("%s\n", cudaGetErrorString (error));
	}
}

void ComputeDensityIISPH_CUDA() {
	ComputeDensityIISPH <<< fcuda.numBlocks, fcuda.numThreads >>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

//void MfComputeDriftVelCUDA ()
//{
//    //mfComputeDriftVel<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
//	
//    error = cudaGetLastError();
//	if (error != cudaSuccess) {
//		fprintf ( stderr, "CUDA ERROR: MfComputeDriftVelCUDA: %s\n", cudaGetErrorString(error) );
//	}    
//	cudaThreadSynchronize ();
//}

//void MfComputeAlphaAdvanceCUDA ()
//{
//	//mfComputeAlphaAdvance<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
//	
//    error = cudaGetLastError();
//	if (error != cudaSuccess) {
//		fprintf ( stderr, "CUDA ERROR: MfComputeAlphaAdvanceCUDA: %s\n", cudaGetErrorString(error) );
//	}    
//	cudaThreadSynchronize ();
//}
//void MfComputeCorrectionCUDA ()
//{
//	//mfComputeCorrection<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );	
//	
//	error = cudaGetLastError();
//	if (error != cudaSuccess) {
//		fprintf ( stderr, "CUDA ERROR: MfComputeCorrectionCUDA: %s\n", cudaGetErrorString(error) );
//	}    
//	cudaThreadSynchronize ();
//}

void MfAdvanceCUDA ()
{
    AdvanceParticles<<< fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum );	
	
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR: MfAdvanceCUDA: %s\n", cudaGetErrorString(error) );
	}    
	cudaThreadSynchronize ();
}

void ComputeForceCUDA_ProjectU(){

	//pressure force, diffusion force
	//ComputeForce_projectu<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	
	error = cudaGetLastError();
	if (error != cudaSuccess)
		fprintf ( stderr, "CUDA ERROR: MfComputeForceCUDA: %s\n", cudaGetErrorString(error) );
	cudaThreadSynchronize ();

	//ComputeSPHtensor<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	
	cudaGetLastError();
	if (error != cudaSuccess)
		fprintf ( stderr, "CUDA ERROR: MfComputSPHtensor: %s\n", cudaGetErrorString(error) );
	cudaThreadSynchronize ();

	//AddSPHtensorForce<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	
	cudaGetLastError();
	if (error != cudaSuccess)
		fprintf ( stderr, "CUDA ERROR: Adding SPH tensor Force: %s\n", cudaGetErrorString(error) );
	cudaThreadSynchronize ();
}

void MfComputeForceCUDA(){

    //SurfaceDetection <<< fcuda.numBlocks, fcuda.numThreads >>> (fbuf, fcuda.pnum);
    //cudaDeviceSynchronize();

	ComputeForce<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	
	error = cudaGetLastError();
	if (error != cudaSuccess)
		fprintf ( stderr, "CUDA ERROR: MfComputeForceCUDA: %s\n", cudaGetErrorString(error) );
	cudaDeviceSynchronize ();

    //SurfaceTension<<< fcuda.numBlocks, fcuda.numThreads >>> (fbuf, fcuda.pnum);
    //error = cudaGetLastError();
	//if (error != cudaSuccess)
	//	fprintf ( stderr, "CUDA ERROR: MfComputeForceCUDA: %s\n", cudaGetErrorString(error) );
    //cudaDeviceSynchronize ();
}



void ComputeMpmForce(){
    
    //Get Grid Mass and Velocity - 1
    //GetGridMassVel <<< fcuda.mpmBlocks, fcuda.mpmThreads>>>(fbuf, fcuda.mpmSize);
    cudaThreadSynchronize();

    //Update Particle Strain Tensor - 2
    //CalcMpmParticleTensor <<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
    cudaThreadSynchronize();

    //Update Grid Force and Velocity - 3
    //CalcMpmGridForce<<<fcuda.mpmBlocks, fcuda.mpmThreads>>>(fbuf, fcuda.mpmSize);
    cudaThreadSynchronize();

    //Update Particle Position and Velocity - 4
    //UpdateMpmParticlePos<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
    cudaThreadSynchronize();

}







//Newly updated 
void ComputeSolidTensorCUDA(){
	//velocity gradient, Strain, Stress
	ComputeSolidTensor<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void ComputeSolidForceCUDA(){
	ComputeSolidForce<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}


void ComputeSolidTensorX_CUDA() {
	//deformation gradient, Strain, Stress
	ComputeSolidTensor_X<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void ComputeSolidForceX_CUDA() {
	ComputeSolidForce_X<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void InitializeSolid_CUDA(){

	GetParticleIndexCUDA();
	GetGridListCUDA();
	RearrageDataCUDA();

	//calculate invA
	ComputeInvA <<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaThreadSynchronize();
}


//IISPH
void ComputeBoundaryDensity() {
	ComputeBoundaryVolume <<< fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void PredictAdvection() {
	//v_adv, dii
	ComputeDii <<<fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();

	//rho_adv_i, aii
	ComputeAii <<<fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void PressureSolve() {
	int iter = 0;
	float rho_avg = 0;

	while (true) {
		Pressure_DP<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
		cudaDeviceSynchronize();

		Pressure_Iter<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
		cudaDeviceSynchronize();

		cudaMemcpy(fbuf.press_l,fbuf.press_l1,sizeof(float)*fcuda.pnum,cudaMemcpyDeviceToDevice);

		//criterion
		thrust::device_ptr<float> d_ptr = thrust::device_pointer_cast(fbuf.densityResidue);
		float sum = thrust::reduce(d_ptr, d_ptr+fcuda.pnum);
		sum /= fcuda.pnum;
		printf("%d iteration: residue %f\n",iter,sum);
		//break;

		if (sum<0.1 && iter>=2) {
			break;
		}
		iter++;
	}


}

void Integration() {
	IntegrateIISPH<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf,fcuda.pnum);
	cudaDeviceSynchronize();
}