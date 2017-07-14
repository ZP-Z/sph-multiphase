

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <conio.h>
//#include <cutil.h>					// cutil32.lib
//#include <cutil_math.h>				// cutil32.lib

//#include <GL/glut.h>
//#include <cuda_gl_interop.h>

#include "fluid_system_host.cuh"		
#include "fluid_system_kern.cuh"
//#include "radixsort.cu"						// Build in RadixSort
#include "thrust\device_vector.h"	//thrust libs
#include "thrust\sort.h" 


FluidParams		fcuda;
//extern ParamCarrier hostCarrier;
bufList			fbuf;
cudaError_t error;

__device__ uint			gridActive;
__device__ int			flagNumFT;  //for transfer
__device__ int			pNumFT;		//for transfer

#define BLOCK_SIZE 256
#define LOCAL_PMAX		896
#define NUM_CELL		27
#define LAST_CELL		26
#define CENTER_CELL		13

float**			g_scanBlockSums;
int**			g_scanBlockSumsInt;
unsigned int	g_numEltsAllocated = 0;
unsigned int	g_numLevelsAllocated = 0;

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
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mpos ) );	
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mvel ) );	
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mveleval ) );	
	CUDA_SAFE_CALL(cudaFree(fbuf.mforce));
	CUDA_SAFE_CALL(cudaFree(fbuf.mpress));
	CUDA_SAFE_CALL(cudaFree(fbuf.last_mpress));
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mdensity ) );		
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mgcell ) );	
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mgndx ) );	
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mclr ) );	
	CUDA_SAFE_CALL ( cudaFree ( fbuf.misbound ) );	

	//CUDA_SAFE_CALL ( cudaFree ( fbuf.mcluster ) );	

	//multi fluid
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_alpha ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_alpha_pre ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_pressure_modify ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_vel_phrel) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_restdensity ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_restmass ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_visc ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_velxcor ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mf_alphagrad ) );

	CUDA_SAFE_CALL ( cudaFree ( fbuf.msortbuf ) );	

	CUDA_SAFE_CALL ( cudaFree ( fbuf.mgrid ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mgridcnt ) );
	//new sort
	CUDA_SAFE_CALL ( cudaFree ( fbuf.midsort ) );

	CUDA_SAFE_CALL ( cudaFree ( fbuf.mgridoff ) );
	CUDA_SAFE_CALL ( cudaFree ( fbuf.mgridactive ) );
}

void FluidSetupCUDA(ParamCarrier& params){
	fcuda.pnum = params.num;
	fcuda.gridRes = params.gridres;
	fcuda.gridSize = params.gridsize;
	fcuda.gridDelta = params.gridIdfac;
	fcuda.gridMin = params.gridmin;
	fcuda.gridMax = params.gridmax;
	fcuda.gridTotal = params.gridtotal;
	fcuda.gridSrch = params.searchnum;
	fcuda.gridAdjCnt = params.neighbornum;
	fcuda.gridScanMax.x = params.gridres.x - params.searchnum;
	fcuda.gridScanMax.y = params.gridres.y - params.searchnum;
	fcuda.gridScanMax.z = params.gridres.z - params.searchnum;
	//fcuda.chk = chk;
	//fcuda.mf_up=0;

	// Build Adjacency Lookup
	int cell = 0;
	for (int y=0; y < params.searchnum; y++)
		for (int z=0; z < params.searchnum; z++)
			for (int x=0; x < params.searchnum; x++)
				fcuda.gridAdj[cell++]  = (y * fcuda.gridRes.z+ z)*fcuda.gridRes.x +  x;

	/*printf ( "CUDA Adjacency Table\n");
	for (int n=0; n < fcuda.gridAdjCnt; n++ ) {
	printf ( "  ADJ: %d, %d\n", n, fcuda.gridAdj[n] );
	}	*/

	// Compute number of blocks and threads
	computeNumBlocks(fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);			// particles
	computeNumBlocks(fcuda.gridTotal, 384, fcuda.gridBlocks, fcuda.gridThreads);		// grid cell
																						// Allocate particle buffers
	fcuda.szPnts = (fcuda.numBlocks  * fcuda.numThreads);
	/*printf ( "CUDA Allocate: \n" );
	printf ( "  Pnts: %d, t:%dx%d=%d, Size:%d\n", fcuda.pnum, fcuda.numBlocks, fcuda.numThreads, fcuda.numBlocks*fcuda.numThreads, fcuda.szPnts);
	printf ( "  Grid: %d, t:%dx%d=%d, bufGrid:%d, Res: %dx%dx%d\n", fcuda.gridTotal, fcuda.gridBlocks, fcuda.gridThreads, fcuda.gridBlocks*fcuda.gridThreads, fcuda.szGrid, (int) fcuda.gridRes.x, (int) fcuda.gridRes.y, (int) fcuda.gridRes.z );
	*/

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mpos, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mvel, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mveleval, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mforce, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)* 3));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mpress, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	//CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.last_mpress, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mdensity, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgcell, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgndx, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mclr, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mColor,    EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float4)));

	int temp_size = EMIT_BUF_RATIO*(4*(sizeof(float)*3) + 2*sizeof(float) + 3*sizeof(uint));

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.misbound, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int)));
	temp_size += EMIT_BUF_RATIO*sizeof(int);
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.accel, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(cfloat3)));
	temp_size += EMIT_BUF_RATIO*sizeof(cfloat3);

	//multi fluid
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_alpha, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*MAX_FLUIDNUM));    //float* num
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_alpha_pre, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*MAX_FLUIDNUM));    //float* num

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_vel_phrel, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3*MAX_FLUIDNUM));	//float*3*num
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_alphagrad, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3*MAX_FLUIDNUM));   //float*3*num

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_pressure_modify, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));				//float
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_restmass, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_restdensity, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_visc, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_velxcor, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3));

	temp_size += EMIT_BUF_RATIO*(2*MAX_FLUIDNUM*sizeof(float) + 2*MAX_FLUIDNUM*(sizeof(float)*3) + 4*sizeof(float) + sizeof(float)*3);

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFtype, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int))); //indicator function
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFid, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int))); //born id
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFidTable, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int))); //id table no sorting

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFtensor, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9)); //deformable tensor
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFRtensor, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFtemptensor, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.MFvelgrad, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9)); //no sorting
	CUDA_SAFE_CALL(cudaMalloc(&fbuf.MFpepsilon, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float))); //no sorting

	temp_size += EMIT_BUF_RATIO*(sizeof(int)*2 + sizeof(float)*9*3);

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mf_multiFlag, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));

	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.msortbuf, EMIT_BUF_RATIO*fcuda.szPnts*temp_size));

	// Allocate grid
	fcuda.szGrid = (fcuda.gridBlocks * fcuda.gridThreads);
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgrid, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridcnt, fcuda.szGrid*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.midsort, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridoff, fcuda.szGrid*sizeof(int)));
	CUDA_SAFE_CALL(cudaMalloc((void**)&fbuf.mgridactive, fcuda.szGrid*sizeof(int)));

	//MpmAllocateBuffer();

	updateParam(&fcuda);
}

void FluidSetupCUDA ( int num, int gsrch, int3 res, cfloat3 size, cfloat3 delta, cfloat3 gmin, cfloat3 gmax, int total, int chk)
{	
	//fcuda.pnum = num;	
	//fcuda.gridRes = res;
	//fcuda.gridSize = size;
	//fcuda.gridDelta = delta;
	//fcuda.gridMin = gmin;
	//fcuda.gridMax = gmax;
	//fcuda.gridTotal = total;
	//fcuda.gridSrch = gsrch;
	//fcuda.gridAdjCnt = gsrch*gsrch*gsrch;
	////fcuda.gridScanMax = res;
	////fcuda.gridScanMax -= make_int3( fcuda.gridSrch, fcuda.gridSrch, fcuda.gridSrch );
	//fcuda.chk = chk;
	//fcuda.mf_up=0;

	//// Build Adjacency Lookup
	//int cell = 0;
	//for (int y=0; y < gsrch; y++ ) 
	//	for (int z=0; z < gsrch; z++ ) 
	//		for (int x=0; x < gsrch; x++ ) 
	//			fcuda.gridAdj [ cell++]  = ( y * fcuda.gridRes.z+ z )*fcuda.gridRes.x +  x ;			
	//
	///*printf ( "CUDA Adjacency Table\n");
	//for (int n=0; n < fcuda.gridAdjCnt; n++ ) {
	//	printf ( "  ADJ: %d, %d\n", n, fcuda.gridAdj[n] );
	//}	*/

	//// Compute number of blocks and threads
	//computeNumBlocks ( fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);			// particles
	//computeNumBlocks ( fcuda.gridTotal, 384, fcuda.gridBlocks, fcuda.gridThreads);		// grid cell
	//// Allocate particle buffers
	//fcuda.szPnts = (fcuda.numBlocks  * fcuda.numThreads);     
	///*printf ( "CUDA Allocate: \n" );
	//printf ( "  Pnts: %d, t:%dx%d=%d, Size:%d\n", fcuda.pnum, fcuda.numBlocks, fcuda.numThreads, fcuda.numBlocks*fcuda.numThreads, fcuda.szPnts);
	//printf ( "  Grid: %d, t:%dx%d=%d, bufGrid:%d, Res: %dx%dx%d\n", fcuda.gridTotal, fcuda.gridBlocks, fcuda.gridThreads, fcuda.gridBlocks*fcuda.gridThreads, fcuda.szGrid, (int) fcuda.gridRes.x, (int) fcuda.gridRes.y, (int) fcuda.gridRes.z );		
	//*/

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mpos,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3 ) );	
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mvel,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3 ) );	
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mveleval,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3 ) );	
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mforce,    EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)* 3));
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mpress,    EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.last_mpress, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)));
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mdensity,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float) ) );	
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mgcell,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mgndx,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint)) );	
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mclr,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint) ) );	
	////CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mColor,    EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float4)));
	//
	//int temp_size = EMIT_BUF_RATIO*(4*(sizeof(float)*3) + 3*sizeof(float) + 2*sizeof(int) + sizeof(uint));

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.misbound, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int)) );	
	//temp_size += EMIT_BUF_RATIO*sizeof(int);
	//CUDA_SAFE_CALL ( cudaMalloc ((void**)&fbuf.accel, EMIT_BUF_RATIO*fcuda.szPnts*sizeof(cfloat3)));
	//temp_size += EMIT_BUF_RATIO*sizeof(cfloat3);

	////multi fluid
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_alpha,					EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*MAX_FLUIDNUM ) );    //float* num
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_alpha_pre,				EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*MAX_FLUIDNUM ) );    //float* num

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_vel_phrel,				EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3*MAX_FLUIDNUM ) );	//float*3*num
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_alphagrad,				EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3*MAX_FLUIDNUM ) );   //float*3*num

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_pressure_modify,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float) ) );				//float
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_restmass,				EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_restdensity,			EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_visc,					EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_velxcor,				EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*3 ) );

	//temp_size += EMIT_BUF_RATIO*(2*MAX_FLUIDNUM*sizeof(float) + 2*MAX_FLUIDNUM*(sizeof(float)*3) + 4*sizeof(float) + sizeof(float)*3 );

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFtype,					EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int) ) ); //indicator function
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFid,						EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int) ) ); //born id
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFidTable,					EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int) ) ); //id table no sorting

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFtensor,			EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9 ) ); //deformable tensor
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFRtensor,			EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9)  );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFtemptensor,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.MFvelgrad,			EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)*9) ); //no sorting
	//CUDA_SAFE_CALL ( cudaMalloc ( &fbuf.MFpepsilon,					EMIT_BUF_RATIO*fcuda.szPnts*sizeof(float)) ); //no sorting

	//temp_size += EMIT_BUF_RATIO*(sizeof(int)*2 + sizeof(float)*9*3 );

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mf_multiFlag,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint) ) );

	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.msortbuf,	EMIT_BUF_RATIO*fcuda.szPnts*temp_size ) );	

	//// Allocate grid
	//fcuda.szGrid = (fcuda.gridBlocks * fcuda.gridThreads);  
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mgrid,		EMIT_BUF_RATIO*fcuda.szPnts*sizeof(int) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mgridcnt,	fcuda.szGrid*sizeof(int) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.midsort,	EMIT_BUF_RATIO*fcuda.szPnts*sizeof(uint) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mgridoff,	fcuda.szGrid*sizeof(int) ) );
	//CUDA_SAFE_CALL ( cudaMalloc ( (void**) &fbuf.mgridactive, fcuda.szGrid*sizeof(int) ) );
	//
 //   //MpmAllocateBuffer();

 //   updateParam( &fcuda );
}

void FluidParamCUDA ( float ss, float sr, float pr, float mass, float rest, cfloat3 bmin, cfloat3 bmax, float estiff, float istiff,float pbstiff, float visc, float damp, float fmin, float fmax, float ffreq, float gslope, float gx, float gy, float gz, float al, float vl )
{
	fcuda.psimscale = ss;
	fcuda.psmoothradius = sr;
	fcuda.pradius = pr;
	fcuda.r2 = sr * sr;
	fcuda.pmass = mass;
	fcuda.prest_dens = rest;	
	fcuda.pboundmin = bmin;
	fcuda.pboundmax = bmax;
	fcuda.pextstiff = estiff;
	fcuda.pintstiff = istiff;
	fcuda.pbstiff = pbstiff;
	fcuda.pvisc = visc;
	fcuda.pdamp = damp;
	fcuda.pforce_min = fmin;
	fcuda.pforce_max = fmax;
	fcuda.pforce_freq = ffreq;
	fcuda.pground_slope = gslope;
	fcuda.pgravity = cfloat3( gx, gy, gz );
	fcuda.AL = al;
	fcuda.AL2 = al * al;
	fcuda.VL = vl;
	fcuda.VL2 = vl * vl;

	printf ( "Bound Min: %f %f %f\n", bmin.x, bmin.y, bmin.z );
	printf ( "Bound Max: %f %f %f\n", bmax.x, bmax.y, bmax.z );

	fcuda.pdist = pow ( fcuda.pmass / fcuda.prest_dens, 1/3.0f );
	fcuda.poly6kern = 315.0f / (64.0f * 3.141592 * pow( sr, 9.0f) );
	fcuda.spikykern = -45.0f / (3.141592 * pow( sr, 6.0f) );
	fcuda.spikykernel = 15 / (3.141592 * pow(sr,6.0f));
	fcuda.lapkern = 45.0f / (3.141592 * pow( sr, 6.0f) );	

	updateParam( &fcuda );
	cudaThreadSynchronize ();
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
	cudaThreadSynchronize();
}

void FluidParamCUDAbuffer_projectu(float* buffer){
	fcuda.coK =					buffer[0];//coK;
	fcuda.coG =					buffer[1];//coG;
	fcuda.phi =					buffer[2];//phi;
	fcuda.coA =					buffer[3];//coA;
	fcuda.coB =					buffer[4];//coB;
	fcuda.coLambdaK =			buffer[5];//coLambdaK;
	fcuda.cohesion =			buffer[6];//cohesion;
	fcuda.boundaryVisc =  		buffer[7];//boundaryVisc;
	fcuda.sleepvel =  			buffer[8];//sleepvel;
	fcuda.initspacing =			buffer[9];//initspacing;
	
	fcuda.coN =					buffer[10];//coN;
	fcuda.Yradius =  			buffer[11];//Yradius;
	fcuda.visc_factor = 		buffer[12];// visc_factor;
	fcuda.fluid_pfactor = 		buffer[13];// fluid_pfactor;
	fcuda.solid_pfactor =  		buffer[14];//solid_pfactor;
	fcuda.fsa = 				buffer[15];// fsa;
	fcuda.fsb =  				buffer[16];//fsb;
	fcuda.bdamp = 				buffer[17];// bdamp;
	fcuda.coD =  				buffer[18];//coD;
	fcuda.coD0 = 				buffer[19];//coD0;
	
	fcuda.solid_coG =  			buffer[20];//solid_coG;
//	fcuda.solid_coV =  			buffer[21];//solid_coV;
	fcuda.solid_coK =  			buffer[22];//solid_coK;
	fcuda.solid_coA =  			buffer[23];//solid_coA;
	fcuda.solid_coB =  			buffer[24];//solid_coB;
	fcuda.solid_fsa =  			buffer[25];//solid_fsa;
	fcuda.solid_fsb =  			buffer[26];//solid_fsb;
	fcuda.solid_coN =  			buffer[27];//solid_coN;
	fcuda.solid_phi =  			buffer[28];//solid_phi;
	fcuda.solid_Yradius =  		buffer[29];//solid_Yradius;
	fcuda.fluidVConstraint =  	buffer[30];//fluidVConstraint;
	fcuda.tohydro =  			buffer[31];//tohydro;

    fcuda.mpmSpacing = buffer[32]; //mpmSpacing
    fcuda.minVec = cfloat3(buffer[33], buffer[34], buffer[35]);
    fcuda.maxVec = cfloat3(buffer[36], buffer[37], buffer[38]);
    
}

void FluidMfParamCUDA ( float *dens, float *visc, float diffusion, float catnum, float dt,  cfloat3 cont, cfloat3 mb1,cfloat3 mb2, float relax,int example)
{
	fcuda.mf_catnum = catnum;
	fcuda.mf_diffusion = diffusion;
	fcuda.mf_dt = dt;
	for(int i=0;i<MAX_FLUIDNUM;i++)
	{
		fcuda.mf_dens[i] = dens[i];
		fcuda.mf_visc[i] = visc[i];
	}
	fcuda.mf_multiFlagPNum = 0;
	//fcuda.mf_splitVolume = splitV;
	//fcuda.mf_mergeVolume = mergeV;
	fcuda.mf_maxPnum = fcuda.pnum * EMIT_BUF_RATIO;
	fcuda.cont =  cont.x;	fcuda.cont1 = cont.y;	fcuda.cont2 = cont.z;	
	fcuda.mb1.x = mb1.x;	fcuda.mb1.y = mb1.y;	fcuda.mb1.z = mb1.z;
	fcuda.mb2.x = mb2.x;	fcuda.mb2.y = mb2.y;	fcuda.mb2.z = mb2.z;
	fcuda.bxmin = mb1.x;    fcuda.by = mb1.y;       fcuda.bzmin = mb1.z;
	fcuda.bxmax = mb2.x;							fcuda.bzmax = mb2.z; 
	
	fcuda.relax = relax;
	fcuda.example = example;
	updateParam( &fcuda );
	cudaThreadSynchronize ();
}


//Copy buffers
void CopyToCUDA ( float* pos, float* vel, float* veleval, float* force, float* pressure, float* density, uint* cluster, uint* gnext, char* clr)
{

	// Send particle buffers
	int numPoints = fcuda.pnum;
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mpos,		pos,			numPoints*sizeof(float)*3, cudaMemcpyHostToDevice ) );	
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mvel,		vel,			numPoints*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mveleval, veleval,		numPoints*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mforce,	force,			numPoints*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mpress,	pressure,		numPoints*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mdensity, density,		numPoints*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mclr,		clr,			numPoints*sizeof(uint), cudaMemcpyHostToDevice ) );


	cudaThreadSynchronize ();	
}
void CopyMfToCUDA ( float* alpha, float* alpha_pre, float* pressure_modify, float* vel_phrel, float* restmass, float* restdensity, float* visc, float* velxcor, float* alphagrad)
{
	// Send particle buffers
	int numPoints = fcuda.pnum;
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_alpha,				alpha,				numPoints*MAX_FLUIDNUM*sizeof(float), cudaMemcpyHostToDevice ) );	
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_alpha_pre,			alpha_pre,			numPoints*MAX_FLUIDNUM*sizeof(float), cudaMemcpyHostToDevice ) );

	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_vel_phrel,			vel_phrel,			numPoints*MAX_FLUIDNUM*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_alphagrad,			alphagrad,			numPoints*MAX_FLUIDNUM*sizeof(float)*3, cudaMemcpyHostToDevice ) );

	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_pressure_modify,	pressure_modify,	numPoints*sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_restmass,			restmass,			numPoints*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_restdensity,		restdensity,		numPoints*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_visc,				visc,				numPoints*sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_velxcor,			velxcor,			numPoints*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	//CUDA_SAFE_CALL( cudaMemcpy ( fbuf.MFtype,				mftype,				numPoints*sizeof(int), cudaMemcpyHostToDevice ) );
	cudaThreadSynchronize ();	
}
void CopyBoundToCUDA (int* isbound )
{
	int numPoints = fcuda.pnum;
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.misbound,	isbound,		numPoints*sizeof(int), cudaMemcpyHostToDevice ) );
	cudaThreadSynchronize ();	
}
void CopyToCUDA_Uproject(int* mftype, float* tensorbuffer, int* bornid){
	int numPoints = fcuda.pnum;
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFtype, mftype, numPoints*sizeof(int), cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFtensor, tensorbuffer, numPoints*sizeof(float)*9, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFtemptensor, tensorbuffer, numPoints*sizeof(float)*9, cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFid, bornid, numPoints*sizeof(int), cudaMemcpyHostToDevice));
	cudaThreadSynchronize ();
}

void CopyFromCUDA ( float* pos, float* vel, float* veleval, float* force, float* pressure, float* density, uint* cluster, uint* gnext, char* clr, int mode)
{
	// Return particle buffers
	int numPoints = fcuda.pnum;

	CUDA_SAFE_CALL( cudaMemcpy ( pos,		fbuf.mpos,			numPoints*sizeof(float)*3, cudaMemcpyDeviceToHost ) );
	CUDA_SAFE_CALL( cudaMemcpy ( clr,		fbuf.mclr,			numPoints*sizeof(uint),  cudaMemcpyDeviceToHost ) );
	CUDA_SAFE_CALL(cudaMemcpy(pressure, fbuf.mpress, numPoints*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL(cudaMemcpy(veleval, fbuf.mforce, numPoints*sizeof(float)*3, cudaMemcpyDeviceToHost));

	if( mode == 2){
		CUDA_SAFE_CALL( cudaMemcpy ( vel,		fbuf.mvel,			numPoints*sizeof(float)*3, cudaMemcpyDeviceToHost ) );
		CUDA_SAFE_CALL( cudaMemcpy ( veleval,	fbuf.mveleval,		numPoints*sizeof(float)*3, cudaMemcpyDeviceToHost ) );
		CUDA_SAFE_CALL( cudaMemcpy ( force,		fbuf.mforce,		numPoints*sizeof(float)*3, cudaMemcpyDeviceToHost ) );
		
		CUDA_SAFE_CALL( cudaMemcpy ( density,	fbuf.mdensity,		numPoints*sizeof(float),  cudaMemcpyDeviceToHost ) );
	}
	
	cudaThreadSynchronize ();	
}
void CopyMfFromCUDA ( float* alpha, float* alpha_pre, float* pressure_modify, float* vel_phrel, float* restmass, float* restdensity, float* visc, float* velxcor, float* alphagrad, int mode)
{

	int numPoints = fcuda.pnum;
	CUDA_SAFE_CALL( cudaMemcpy ( alpha,				fbuf.mf_alpha,				numPoints*MAX_FLUIDNUM*sizeof(float), cudaMemcpyDeviceToHost ) );
	CUDA_SAFE_CALL( cudaMemcpy ( restmass,			fbuf.mf_restmass,			numPoints*sizeof(float),  cudaMemcpyDeviceToHost ) );
	CUDA_SAFE_CALL( cudaMemcpy ( restdensity,		fbuf.mf_restdensity,		numPoints*sizeof(float),  cudaMemcpyDeviceToHost ) );
	if( mode == 2){
		// Send particle buffers
		
		CUDA_SAFE_CALL( cudaMemcpy ( alpha_pre,			fbuf.mf_alpha_pre,			numPoints*MAX_FLUIDNUM*sizeof(float), cudaMemcpyDeviceToHost ) );
		CUDA_SAFE_CALL( cudaMemcpy ( pressure_modify,	fbuf.mf_pressure_modify,	numPoints*sizeof(float), cudaMemcpyDeviceToHost ) );
		CUDA_SAFE_CALL( cudaMemcpy ( vel_phrel,			fbuf.mf_vel_phrel,			numPoints*MAX_FLUIDNUM*sizeof(float)*3, cudaMemcpyDeviceToHost ) );

		CUDA_SAFE_CALL( cudaMemcpy ( visc,				fbuf.mf_visc,				numPoints*sizeof(float), cudaMemcpyDeviceToHost ) );
		CUDA_SAFE_CALL( cudaMemcpy ( velxcor,			fbuf.mf_velxcor,			numPoints*sizeof(float)*3, cudaMemcpyDeviceToHost ) );
		CUDA_SAFE_CALL( cudaMemcpy ( alphagrad,			fbuf.mf_alphagrad,			numPoints*MAX_FLUIDNUM*sizeof(float)*3, cudaMemcpyDeviceToHost ) );
	}
}
void CopyBoundFromCUDA (int* isbound )
{
	int numPoints = fcuda.pnum;
	if ( isbound != 0x0 ) CUDA_SAFE_CALL( cudaMemcpy ( isbound,	fbuf.misbound,		numPoints*sizeof(int),  cudaMemcpyDeviceToHost ) );
	cudaThreadSynchronize ();	
}
void CopyFromCUDA_Uproject(int* mftype, int* idtable, float* pepsilon, float* stensor, int mode){
	int numPoints = fcuda.pnum;
	CUDA_SAFE_CALL( cudaMemcpy( mftype, fbuf.MFtype, numPoints*sizeof(int), cudaMemcpyDeviceToHost));
	CUDA_SAFE_CALL( cudaMemcpy( idtable, fbuf.MFid, numPoints*sizeof(int), cudaMemcpyDeviceToHost));
	//CUDA_SAFE_CALL( cudaMemcpy( pepsilon, fbuf.MFpepsilon, numPoints*sizeof(float), cudaMemcpyDeviceToHost));
	
	if( mode==2) //to save particles
		CUDA_SAFE_CALL( cudaMemcpy( stensor, fbuf.MFtemptensor, numPoints*sizeof(float)*9, cudaMemcpyDeviceToHost));
	cudaThreadSynchronize ();
}


//Called when particles emitted
void CopyEmitToCUDA ( float* pos, float* vel, float* veleval, float* force, float* pressure, float* density, uint* cluster, uint* gnext, char* clr, int startnum, int numcount,int* isbound )
{

	// Send particle buffers
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mpos+startnum,		pos+startnum*3,			numcount*sizeof(float)*3, cudaMemcpyHostToDevice ) );	
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mvel+startnum,		vel+startnum*3,			numcount*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mveleval+startnum,	veleval+startnum*3,		numcount*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mforce+startnum,	force+startnum*3,			numcount*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mpress+startnum,		pressure+startnum,		numcount*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mdensity+startnum,	density+startnum,		numcount*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mclr+startnum,		clr+startnum,			numcount*sizeof(uint), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.misbound + startnum,				isbound + startnum,				numcount*sizeof(int), cudaMemcpyHostToDevice ) );

	cudaThreadSynchronize ();	
}	
void CopyEmitMfToCUDA ( float* alpha, float* alpha_pre, float* pressure_modify, float* vel_phrel, float* restmass, float* restdensity, float* visc, float* velxcor, float* alphagrad,int startnum, int numcount)
{
	// Send particle buffers
	int mulstartnum = startnum*MAX_FLUIDNUM;
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_alpha + mulstartnum,				alpha + mulstartnum,				numcount*MAX_FLUIDNUM*sizeof(float), cudaMemcpyHostToDevice ) );	
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_alpha_pre + mulstartnum,			alpha_pre + mulstartnum,			numcount*MAX_FLUIDNUM*sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_pressure_modify+startnum,			pressure_modify+startnum,			numcount*sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_vel_phrel + mulstartnum,			vel_phrel + mulstartnum*3,			numcount*MAX_FLUIDNUM*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_restmass+startnum,					restmass+startnum,					numcount*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_restdensity+startnum,				restdensity+startnum,				numcount*sizeof(float),  cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_visc+startnum,						visc+startnum,						numcount*sizeof(float), cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_velxcor+startnum,					velxcor+startnum*3,					numcount*sizeof(float)*3, cudaMemcpyHostToDevice ) );
	CUDA_SAFE_CALL( cudaMemcpy ( fbuf.mf_alphagrad + mulstartnum,			alphagrad + mulstartnum*3,			numcount*MAX_FLUIDNUM*sizeof(float)*3, cudaMemcpyHostToDevice ) );

	cudaThreadSynchronize ();	
}

void CopyEmitToCUDA_Uproject(int* mftype, float* tensorbuffer, int* bornid, int startnum, int numcount){
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFtype + startnum,			mftype+startnum,				numcount*sizeof(int),			cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFtensor + startnum*9,		tensorbuffer+startnum*9,		numcount*sizeof(float)*9,		cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFtemptensor + startnum*9, tensorbuffer+startnum*9,		numcount*sizeof(float)*9,		cudaMemcpyHostToDevice));
	CUDA_SAFE_CALL( cudaMemcpy( fbuf.MFid + startnum,			bornid+startnum,				numcount*sizeof(int),			cudaMemcpyHostToDevice));
	
	cudaThreadSynchronize ();
}

void UpdatePNumCUDA( int newPnum)
{
	fcuda.pnum = newPnum;
	computeNumBlocks ( fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);    //threads changed!
	fcuda.szPnts = (fcuda.numBlocks  * fcuda.numThreads);					   //szPnts changed!	
	updateParam( &fcuda );
	cudaThreadSynchronize ();
}
int MfGetPnum(){
	return fcuda.pnum;
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
	cudaMemcpy ( fbuf.msortbuf + n*BUF_POS,		fbuf.mpos,		n*sizeof(float)*3,	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_VEL,		fbuf.mvel,		n*sizeof(float)*3,	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_VELEVAL,	fbuf.mveleval,	n*sizeof(float)*3,	cudaMemcpyDeviceToDevice );
	cudaMemcpy(fbuf.msortbuf + n*BUF_FORCE,		fbuf.mforce, n*sizeof(float)* 3, cudaMemcpyDeviceToDevice);
	cudaMemcpy(fbuf.msortbuf + n*BUF_PRESS,		fbuf.mpress, n*sizeof(float), cudaMemcpyDeviceToDevice);
	//cudaMemcpy(fbuf.msortbuf + n*BUF_LAST_PRESS, fbuf.last_mpress, n*sizeof(float), cudaMemcpyDeviceToDevice);
	cudaMemcpy ( fbuf.msortbuf + n*BUF_DENS,	fbuf.mdensity,	n*sizeof(float),	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_GCELL,	fbuf.mgcell,	n*sizeof(uint),		cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_GNDX,	fbuf.mgndx,		n*sizeof(uint),		cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_CLR,		fbuf.mclr,		n*sizeof(uint),		cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_ISBOUND,		fbuf.misbound,		n*sizeof(int),		cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_ACCEL,	fbuf.accel,		n*sizeof(cfloat3),	cudaMemcpyDeviceToDevice);

	//multi fluid
	cudaMemcpy ( fbuf.msortbuf + n*BUF_ALPHA,			fbuf.mf_alpha,			n*MAX_FLUIDNUM*sizeof(float),	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_ALPHAPRE,		fbuf.mf_alpha_pre,		n*MAX_FLUIDNUM*sizeof(float),	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_PRESSMODI,		fbuf.mf_pressure_modify,n*sizeof(float),				cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_VELPHREL,		fbuf.mf_vel_phrel,		n*MAX_FLUIDNUM*sizeof(cfloat3),	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_RMASS,			fbuf.mf_restmass,		n*sizeof(float),				cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_RDENS,			fbuf.mf_restdensity,	n*sizeof(float),				cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_VISC,			fbuf.mf_visc,			n*sizeof(float),				cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_VELXCOR,			fbuf.mf_velxcor,		n*sizeof(cfloat3),				cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_ALPHAGRAD,		fbuf.mf_alphagrad,		n*MAX_FLUIDNUM*sizeof(cfloat3),	cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_INDICATOR,		fbuf.MFtype,			n*sizeof(int),					cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_TENSOR,			fbuf.MFtensor,		    n*sizeof(float)*9,				cudaMemcpyDeviceToDevice );	
	cudaMemcpy ( fbuf.msortbuf + n*BUF_TEMPTENSOR,      fbuf.MFtemptensor,		n*sizeof(float)*9,				cudaMemcpyDeviceToDevice );	
	cudaMemcpy ( fbuf.msortbuf + n*BUF_RTENSOR,			fbuf.MFRtensor,			n*sizeof(float)*9,				cudaMemcpyDeviceToDevice );
	cudaMemcpy ( fbuf.msortbuf + n*BUF_BORNID,			fbuf.MFid,				n*sizeof(int),					cudaMemcpyDeviceToDevice );

	// Counting Sort - pass one, determine grid counts
	cudaMemset ( fbuf.mgrid,	GRID_UCHAR,	fcuda.pnum * sizeof(int) );
	
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
	////CUDA_SAFE_CALL(cudaMemcpyFromSymbol(&(fcuda.pnum), pNumFT, sizeof(int)));  //total pnum changed!
	////computeNumBlocks ( fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);    //threads changed!
	////fcuda.szPnts = (fcuda.numBlocks  * fcuda.numThreads);					   //szPnts changed!
	////		printf("pnum:%d,Blocknum:%d,Threadnum:%d\n",fcuda.pnum,fcuda.numBlocks,fcuda.numThreads);
	////cudaThreadSynchronize ();
}

void MpmAllocateBuffer(){

    //calculate size
    int xlen = (fcuda.maxVec - fcuda.minVec).x / fcuda.mpmSpacing;
    int ylen = (fcuda.maxVec - fcuda.minVec).y / fcuda.mpmSpacing;
    int zlen = (fcuda.maxVec - fcuda.minVec).z / fcuda.mpmSpacing;
    fcuda.mpmSize = xlen*ylen*zlen;
    fcuda.mpmXl = xlen;
    fcuda.mpmYl = ylen;
    fcuda.mpmZl = zlen;

    int splitnum = 2;

    computeNumBlocks ( fcuda.mpmSize, 384, fcuda.mpmBlocks, fcuda.mpmThreads);

    cudaMalloc(&fbuf.mpmMass,   fcuda.mpmSize * sizeof(float) * splitnum);
    cudaMalloc(&fbuf.mpmVel,    fcuda.mpmSize * sizeof(cfloat3) * splitnum); //0-solid, 1-fluid
    
    cudaMalloc(&fbuf.mpmAlpha, fcuda.mpmSize * sizeof(float) * MAX_FLUIDNUM);
    cudaMalloc(&fbuf.mpmForce, fcuda.mpmSize * sizeof(cfloat3) * splitnum);

    cudaMalloc(&fbuf.mpmTensor, fcuda.mpmSize * sizeof(float) * 9);
    cudaMalloc(&fbuf.mpmGid,    fcuda.mpmSize * sizeof(uint));
    cudaMalloc(&fbuf.mpmGridVList, fcuda.mpmSize * sizeof(uint));
    cudaMalloc(&fbuf.mpmIdSort, fcuda.mpmSize * sizeof(uint));
    cudaMalloc(&fbuf.mpmGridCnt, fcuda.szGrid * sizeof(uint));
    cudaMalloc(&fbuf.mpmGridOff, fcuda.szGrid * sizeof(uint));

    cudaMalloc(&fbuf.mpmPos,    fcuda.mpmSize * sizeof(cfloat3));
}

void MpmSortGridCuda(){
    //initMpm <<< fcuda.mpmBlocks, fcuda.mpmThreads >>> (fbuf, fcuda.mpmSize);
    cudaThreadSynchronize();

	thrust::device_ptr<uint> dev_keysg ( fbuf.mpmGid );
	thrust::device_ptr<uint> dev_valuesg ( fbuf.mpmIdSort );
	thrust::sort_by_key ( dev_keysg, dev_keysg + fcuda.mpmSize, dev_valuesg );
	cudaThreadSynchronize ();
	
    //MpmCalcFirstCnt << < fcuda.mpmBlocks, fcuda.mpmThreads >> > (fbuf, fcuda.mpmSize);
	cudaThreadSynchronize ();

	//MpmGetCnt << <fcuda.mpmBlocks, fcuda.mpmThreads >> > (fbuf, fcuda.mpmSize);
	cudaThreadSynchronize ();
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
	
	//ComputeBoundaryVolume <<< fcuda.numBlocks, fcuda.numThreads>>> (fbuf, fcuda.pnum);
	//cudaDeviceSynchronize();

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
	cudaMemset(fbuf.MFtensor, 0, sizeof(float)*9*fcuda.pnum);
	cudaMemset(fbuf.accel, 0, sizeof(cfloat3)*fcuda.pnum);
}