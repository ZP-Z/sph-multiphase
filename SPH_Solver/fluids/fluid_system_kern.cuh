
#ifndef DEF_KERN_CUDA
	#define DEF_KERN_CUDA

#include "..\\fluids\\multifluid_def.h"

#include <stdio.h>
#include <math.h>

#include "../common_header.h"
#include "geometry.h"
	

	// Temporary sort buffer offsets
	
#define BUF_DISPLAYBUF 0
#define BUF_CALCBUF (BUF_DISPLAYBUF + sizeof(displayPack))
#define BUF_INTMBUF (BUF_CALCBUF + sizeof(calculationPack))

	// metadata
	struct FluidParams {
		int	pnum;
		int	szPnts;
		int	numThreads, numBlocks;

		int	gridTotal;
		int szGrid;
		int	gridThreads, gridBlocks;	

        int mpmSize;
        int mpmBlocks, mpmThreads;
    };

	struct bufList {
		//Particle properties

		displayPack* displayBuffer;
		calculationPack* calcBuffer;
		IntermediatePack* intmBuffer;

		//For sorting
		char*			msortbuf;
		int*			mgridcnt;
		int*			mgridoff;
		uint*			mgcell;
		uint*			mgndx;
		int*			MFidTable;
		uint*			midsort;
		//End sorting

		//Mpm/Flip Grid
		cfloat3*        mpmPos; //node
		float*          mpmMass;
		cfloat3*        mpmVel; //node velocity
		float *u, *v, *w; //staggered velocity
		cfloat3*        mpmForce;
		uint*           mpmGid;    //mpmSize

	};// End particle&grid buffers
	
	//new sort
	__global__ void InitialSort ( bufList buf, int pnum );
	__global__ void CalcFirstCnt ( bufList buf, int pnum );
	__global__ void RearrangeData ( bufList buf, int pnum );
	__global__ void GetCnt ( bufList buf, int gnum );

	__global__ void initDensity(bufList buf,int pnum);

	//multi fluid
	//calculating functions

	__global__ void mfFindNearest (bufList buf,int pnum);

	__global__ void mfComputeDriftVel( bufList buf, int pnum );
	__global__ void mfComputeAlphaAdvance( bufList buf, int pnum );
	__global__ void mfComputeCorrection( bufList buf, int pnum );

	__global__ void AdvanceParticles(bufList buf, int numPnts );
    

	__global__ void ComputeDensityPressure(bufList buf,int pnum);
	__global__ void ComputeBoundaryVolume(bufList buf,int pnum);



	//calculating functions for project-u
	__global__ void ComputeForce_projectu( bufList buf, int pnum );

	__global__ void ComputeSPHtensor(bufList buf,int pnum);
	__global__ void AddSPHtensorForce(bufList buf,int pnum);
	//end calculating functions for project-u

    //Surface Tension Force
    __global__ void SurfaceDetection(bufList buf,int pnum);
    __global__ void ComputeForce( bufList buf, int pnum );
    __global__ void SurfaceTension(bufList buf,int pnum);


    //Mpm 

 
    __global__ void GetGridMassVel(bufList buf,int mpmSize);
    
	__global__ void CalcMpmParticleTensor(bufList buf,int pnum);
    __global__ void CalcMpmGridForce(bufList buf,int mpmSize);
    __global__ void UpdateMpmParticlePos(bufList buf, int pnum);
    
	//newly updated solid functions
	
	__global__ void ComputeSolidTensor_CUDA(bufList buf,int pnum);
	__global__ void ComputeSolidForce_CUDA(bufList buf,int pnum);
	__global__ void ComputeDensity_CUDA(bufList buf, int pnum);

	void CarryParam(ParamCarrier& hostCarrier);

	#define EPSILON				0.00001f
	#define GRID_UCHAR			0xFF
	#define GRID_UNDEF			4294967295

#endif
