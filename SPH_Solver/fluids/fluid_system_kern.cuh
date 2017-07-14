
#ifndef DEF_KERN_CUDA
	#define DEF_KERN_CUDA

#include "..\\fluids\\multifluid_def.h"

#include <stdio.h>
#include <math.h>

#include "../common_header.h"
#include "../geometry.h"
	

	// Temporary sort buffer offsets
	
#define BUF_DISPLAYBUF 0
#define BUF_CALCBUF (BUF_DISPLAYBUF + sizeof(displayPack))
#define BUF_INTMBUF (BUF_CALCBUF + sizeof(calculationPack))

	// Fluid Parameters (stored on both host and device)
	struct FluidParams {
		int				numThreads, numBlocks;
		int				gridThreads, gridBlocks;	

		int				szPnts, szHash, szGrid;
		int				stride, pnum;
		int				chk;
		float			pdist, pmass, prest_dens;
		float			pextstiff, pintstiff, pbstiff;
		float			pradius, psmoothradius, r2, psimscale, pvisc;
		float			pforce_min, pforce_max, pforce_freq, pground_slope;
		float			pvel_limit, paccel_limit, pdamp;
		cfloat3			pboundmin, pboundmax, pgravity; //p is soft bound
		cfloat3			mb1,mb2;
		float			AL, AL2, VL, VL2; //limits of acceleration and velocity
		
		float			poly6kern, spikykern, lapkern;

		float spikykernel;//not the derivative

		cfloat3			gridSize, gridDelta, gridMin, gridMax;
		cint3			gridRes, gridScanMax;
		int				gridSrch, gridTotal, gridAdjCnt, gridActive;
		float			test1,test2,test3;
		int				gridAdj[64];
		float			coLambdaK, cohesion;
		//multi fluid parameters
		float			mf_dens[MAX_FLUIDNUM];
		float			mf_visc[MAX_FLUIDNUM];
		float			mf_diffusion;
		int				mf_catnum;
		float			mf_dt;

		uint			mf_multiFlagPNum;   //used with Buflist.mf_multiFlag, stores total number count (of flags)
		float			mf_splitVolume;
		float			mf_mergeVolume;
		uint			mf_maxPnum;
		float			cont,cont1,cont2;
		int				mf_up;
		float relax;
		int example;
		float			by,bxmin,bxmax,bzmin,bzmax,pan_r,omega; // for example3 rotation
		int				loadwhich;
		float			solid_coG, solid_coV, solid_coK, solid_coA, solid_coB, solid_fsa, solid_fsb, solid_coN, solid_phi, solid_Yradius;

		float			coK, coG, phi, coA, coB, initspacing, coN, Yradius;
		float			visc_factor, solid_pfactor, fluid_pfactor;
		float			fsa, fsb, bdamp, boundaryVisc, sleepvel,coD,coD0;
		float			fluidVConstraint, tohydro;

        int             mpmXl, mpmYl, mpmZl;
        int             mpmSize;
        int             mpmBlocks, mpmThreads;
        float           mpmSpacing;
        cfloat3          minVec, maxVec;

    };

	

	// Prefix Sum defines - 16 banks on G80
	#define NUM_BANKS		16
	#define LOG_NUM_BANKS	 4
	
	//new sort
	__global__ void InitialSort ( bufList buf, int pnum );
	__global__ void CalcFirstCnt ( bufList buf, int pnum );
	__global__ void CountingSortFull_ ( bufList buf, int pnum );
	__global__ void GetCnt ( bufList buf, int gnum );

	__global__ void initDensity(bufList buf,int pnum);

	//multi fluid
	//calculating functions

	__global__ void mfFindNearest (bufList buf,int pnum);

	__global__ void mfComputeDriftVel( bufList buf, int pnum );
	__global__ void mfComputeAlphaAdvance( bufList buf, int pnum );
	__global__ void mfComputeCorrection( bufList buf, int pnum );

	__global__ void AdvanceParticles( float time, float dt, float ss, bufList buf, int numPnts );
    

	__global__ void ComputeDensityPressure(bufList buf,int pnum);
	__global__ void ComputeBoundaryVolume(bufList buf,int pnum);

	//calculating functions for certain cases
	__global__ void mfChangeDensity (bufList buf,int pnum,const float scale);

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

    __global__ void initMpm        (bufList buf, int mpmSize);
    __global__ void MpmGetCnt      (bufList buf,int mpmSize);
    __global__ void MpmCalcFirstCnt(bufList buf,int mpmSize);

    __global__ void GetGridMassVel(bufList buf,int mpmSize);
    
	__global__ void CalcMpmParticleTensor(bufList buf,int pnum);
    __global__ void CalcMpmGridForce(bufList buf,int mpmSize);
    __global__ void UpdateMpmParticlePos(bufList buf, int pnum);
    
	//newly updated solid functions
	
	__global__ void ComputeSolidTensor_CUDA(bufList buf,int pnum);
	__global__ void ComputeSolidForce_CUDA(bufList buf,int pnum);
	__global__ void ComputeDensity_CUDA(bufList buf, int pnum);

    void updateParam( FluidParams* paramCPU );
	void CarryParam(ParamCarrier& hostCarrier);

	#define EPSILON				0.00001f
	#define GRID_UCHAR			0xFF
	#define GRID_UNDEF			4294967295

#endif
