

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <conio.h>


#include "fluid_system_host.cuh"		
#include "fluid_system_kern.cuh"
#include "MpmSolver.cuh"

#include "thrust\device_vector.h"	//thrust libs
#include "thrust\sort.h" 


FluidParams		fcuda;
bufList			fbuf;

cudaError_t error;

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
	cudaFree(fbuf.intmBuffer);

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
	fcuda.gridRes = params.gridres;
	fcuda.gridSize = params.gridsize;
	fcuda.gridDelta = params.gridIdfac;
	fcuda.gridMin = params.gridmin;
	fcuda.gridMax = params.gridmax;
	fcuda.gridTotal = params.gridtotal;
	fcuda.gridSrch = params.searchnum; //3
	fcuda.gridAdjCnt = params.neighbornum;
	fcuda.gridScanMax.x = params.gridres.x - params.searchnum;
	fcuda.gridScanMax.y = params.gridres.y - params.searchnum;
	fcuda.gridScanMax.z = params.gridres.z - params.searchnum;
	//fcuda.chk = chk;
	//fcuda.mf_up=0;

	// Build Adjacency Lookup
	int cell = 0;

	for (int y=-1; y <= 1; y++)
		for (int z=-1; z <=1; z++)
			for (int x=-1; x <= 1; x++)
				fcuda.gridAdj[cell++]  = y*fcuda.gridRes.z*fcuda.gridRes.x + z*fcuda.gridRes.x +  x;

	/*printf ( "CUDA Adjacency Table\n");
	for (int n=0; n < fcuda.gridAdjCnt; n++ ) {
	printf ( "  ADJ: %d, %d\n", n, fcuda.gridAdj[n] );
	}	*/

	// Compute number of blocks and threads
	computeNumBlocks(fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);			// particles
	computeNumBlocks(fcuda.gridTotal, 384, fcuda.gridBlocks, fcuda.gridThreads);		// grid cell
	
    /*printf ( "CUDA Allocate: \n" );
	printf ( "  Pnts: %d, t:%dx%d=%d, Size:%d\n", fcuda.pnum, fcuda.numBlocks, fcuda.numThreads, fcuda.numBlocks*fcuda.numThreads, fcuda.szPnts);
	printf ( "  Grid: %d, t:%dx%d=%d, bufGrid:%d, Res: %dx%dx%d\n", fcuda.gridTotal, fcuda.gridBlocks, fcuda.gridThreads, fcuda.gridBlocks*fcuda.gridThreads, fcuda.szGrid, (int) fcuda.gridRes.x, (int) fcuda.gridRes.y, (int) fcuda.gridRes.z );
	*/


	// Allocate particle buffers
	//fcuda.szPnts = (fcuda.numBlocks  * fcuda.numThreads);
	fcuda.szPnts = params.maxNum;

	cudaMalloc(&fbuf.displayBuffer, EMIT_BUF_RATIO * fcuda.szPnts * sizeof(displayPack));
	cudaMalloc(&fbuf.calcBuffer,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(calculationPack));
	cudaMalloc(&fbuf.intmBuffer,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(IntermediatePack));
	
	int temp_size = EMIT_BUF_RATIO*(sizeof(displayPack) + sizeof(calculationPack));


	cudaMalloc(&fbuf.MFidTable,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int)); //id table no sorting
	cudaMalloc(&fbuf.msortbuf,		EMIT_BUF_RATIO*fcuda.szPnts*temp_size);

	// Allocate grid
	fcuda.szGrid = (fcuda.gridBlocks * fcuda.gridThreads);
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgcell, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgndx, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridcnt, fcuda.szGrid*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.midsort, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridoff, fcuda.szGrid*sizeof(int)));
	

	//MpmAllocateBuffer();

	updateParam(&fcuda);
}


void FluidParamCUDA (ParamCarrier& params){
	fcuda.psimscale = params.simscale;
	fcuda.psmoothradius = params.smoothradius; //real smooth radius
	fcuda.pradius = params.radius;
	fcuda.r2 = params.smoothradius * params.smoothradius;
	fcuda.pmass = params.mass;
	fcuda.prest_dens = params.restdensity;
	fcuda.pvisc = params.viscosity;
	fcuda.pboundmin = params.softminx;
	fcuda.pboundmax = params.softmaxx;
	fcuda.pextstiff = params.extstiff;
	fcuda.pintstiff = params.intstiff;
//	fcuda.pbstiff = pbstiff;
	fcuda.pdamp = params.extdamp;
	//fcuda.pforce_min = fmin;
	//fcuda.pforce_max = fmax;
	//fcuda.pforce_freq = ffreq;
	//fcuda.pground_slope = gslope;
	fcuda.pgravity = params.gravity;
	fcuda.AL = params.acclimit;
	fcuda.AL2 = params.acclimit * params.acclimit;
	fcuda.VL = params.vlimit;
	fcuda.VL2 = params.vlimit * params.vlimit;

	printf("Bound Min: %f %f %f\n", fcuda.pboundmin.x, fcuda.pboundmin.y, fcuda.pboundmin.z);
	printf("Bound Max: %f %f %f\n", fcuda.pboundmax.x, fcuda.pboundmax.y, fcuda.pboundmax.z);

	fcuda.pdist = pow(fcuda.pmass / fcuda.prest_dens, 1/3.0f);
	fcuda.poly6kern = 315.0f / (64.0f * 3.141592 * pow(fcuda.psmoothradius, 9.0f));
	fcuda.spikykern = -45.0f / (3.141592 * pow(fcuda.psmoothradius, 6.0f));
	fcuda.spikykernel = 15 / (3.141592 * pow(fcuda.psmoothradius, 6.0f));
	fcuda.lapkern = 45.0f / (3.141592 * pow(fcuda.psmoothradius, 6.0f));

	//fcuda.mf_catnum = catnum;
	//fcuda.mf_diffusion = diffusion;
	fcuda.mf_dt = params.dt;
	/*for (int i=0; i<MAX_FLUIDNUM; i++)
	{
		fcuda.mf_dens[i] = dens[i];
		fcuda.mf_visc[i] = visc[i];
	}*/

	updateParam(&fcuda);
}






//Called in RunSimulateCudaFull
void InitialSortCUDA( uint* gcell, uint* ccell, int* gcnt )
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

	// Transfer data back if requested (for validation)
	if (gcell != 0x0) {
		CUDA_SAFE_CALL( cudaMemcpy ( gcell,	fbuf.mgcell,	fcuda.pnum*sizeof(uint),		cudaMemcpyDeviceToHost ) );		
		CUDA_SAFE_CALL( cudaMemcpy ( gcnt,	fbuf.mgridcnt,	fcuda.gridTotal*sizeof(int),	cudaMemcpyDeviceToHost ) );
	}
}
void SortGridCUDA( int* goff )
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

void CountingSortFullCUDA_( uint* ggrid )
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

	CountingSortFull_ <<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );		
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
    //
    //cudaMalloc(&fbuf.mpmAlpha, fcuda.mpmSize * sizeof(float) * MAX_FLUIDNUM);
    //cudaMalloc(&fbuf.mpmForce, fcuda.mpmSize * sizeof(cfloat3) * splitnum);

	//cint3 tmp = param.mpmRes;
	//cudaMalloc(&fbuf.u,	(tmp.x+1)*tmp.y*tmp.z*sizeof(float));
	//cudaMalloc(&fbuf.v, (tmp.x+1)*tmp.y*tmp.z*sizeof(float));
	//cudaMalloc(&fbuf.w, (tmp.x+1)*tmp.y*tmp.z*sizeof(float));

    //cudaMalloc(&fbuf.mpmTensor, fcuda.mpmSize * sizeof(float) * 9);
    cudaMalloc(&fbuf.mpmGid,		fcuda.mpmSize * sizeof(uint));

	updateParam(&fcuda);
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
	cudaThreadSynchronize();

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
	
	/*mfFindNearest<<< fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR: MfFindNearestVelCUDA: %s\n", cudaGetErrorString(error) );
	}    
	cudaDeviceSynchronize ();*/
	
	ComputeBoundaryVolume <<< fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();

	ComputeDensityPressure <<< fcuda.numBlocks, fcuda.numThreads >>> (fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
	
	error = cudaGetLastError();
	if(error != cudaSuccess){
		printf("%s\n", cudaGetErrorString (error));
	}
}

void MfComputeDriftVelCUDA ()
{
    //mfComputeDriftVel<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	
    error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR: MfComputeDriftVelCUDA: %s\n", cudaGetErrorString(error) );
	}    
	cudaThreadSynchronize ();
}

void MfComputeAlphaAdvanceCUDA ()
{
	//mfComputeAlphaAdvance<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );
	
    error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR: MfComputeAlphaAdvanceCUDA: %s\n", cudaGetErrorString(error) );
	}    
	cudaThreadSynchronize ();
}
void MfComputeCorrectionCUDA ()
{
	//mfComputeCorrection<<< fcuda.numBlocks, fcuda.numThreads>>> ( fbuf, fcuda.pnum );	
	
	error = cudaGetLastError();
	if (error != cudaSuccess) {
		fprintf ( stderr, "CUDA ERROR: MfComputeCorrectionCUDA: %s\n", cudaGetErrorString(error) );
	}    
	cudaThreadSynchronize ();
}

void MfAdvanceCUDA ( float time , float dt, float ss )
{
    AdvanceParticles<<< fcuda.numBlocks, fcuda.numThreads>>> ( time, dt, ss, fbuf, fcuda.pnum );	
	
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
void ComputeSolidTensor(){
	//Get Density
	//ComputeDensity_CUDA<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();

	//velocity gradient, Strain, Stress
	//ComputeSolidTensor_CUDA<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void ComputeSolidForce(){
	//ComputeSolidForce_CUDA<<<fcuda.numBlocks, fcuda.numThreads>>>(fbuf, fcuda.pnum);
	cudaDeviceSynchronize();
}

void InitSolid(){
	//cudaMemset(fbuf.MFtensor, 0, sizeof(float)*9*fcuda.pnum);
	//cudaMemset(fbuf.accel, 0, sizeof(cfloat3)*fcuda.pnum);
}