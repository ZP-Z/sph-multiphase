/*
  FLUIDS v.3 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2012. Rama Hoetzlein, http://fluids3.com

  Fluids-ZLib license (* see part 1 below)
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
	 claim that you wrote the original software. Acknowledgement of the
	 original author is required if you publish this in a paper, or use it
	 in a product. (See fluids3.com for details)
  2. Altered source versions must be plainly marked as such, and must not be
	 misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
*/
	 
#ifndef DEF_HOST_CUDA
#define DEF_HOST_CUDA

#include "..\\fluids\\multifluid_def.h"
#include <cuda.h>
#include <cuda_runtime_api.h>
#include <cuda_runtime.h>
#include <device_launch_parameters.h>
#include "fluid_system_kern.cuh"

	//#define TOTAL_THREADS			1000000
	//#define BLOCK_THREADS			256
	//#define MAX_NBR					80	
	
	//#define COLOR(r,g,b)	( (uint((b)*255.0f)<<16) | (uint((g)*255.0f)<<8) | uint((r)*255.0f) )
	#define COLORA(r,g,b,a)	( (uint((a)*255.0f)<<24) | (uint((b)*255.0f)<<16) | (uint((g)*255.0f)<<8) | uint((r)*255.0f) )

	typedef unsigned int		uint;
	typedef unsigned short		ushort;
	typedef unsigned char		uchar;

	extern "C"
	{

	void cudaInit(int argc, char **argv);
	void cudaExit(int argc, char **argv);

	void FluidClearCUDA ();


	void FluidSetupCUDA(ParamCarrier& params);
	void FluidParamCUDA(ParamCarrier& param);

	//new sort
	void InitialSortCUDA( uint* gcell, uint* ccell, int* gcnt );
	void SortGridCUDA( int* goff );
	void CountingSortFullCUDA_( uint* ggrid );
	void initSPH(float* restdensity,int* mftype);
	void InitSolid();

	//Multifluid simulation
	void MfComputePressureCUDA();
	void MfComputeDriftVelCUDA();
	void MfComputeAlphaAdvanceCUDA();
	void MfComputeCorrectionCUDA();  
	void MfComputeForceCUDA ();	
	void MfAdvanceCUDA ( float time , float dt, float ss );

	//Project-U changes computing force
	void ComputeForceCUDA_ProjectU();

    //Mpm
	void MpmAllocateBufferCUDA(ParamCarrier& param);
    void ComputeMpmForce();
	void IndexMPMSortCUDA();
	void MpmColorTestCUDA();
	void MpmGetMomentumCUDA();


	void ComputeSolidTensorCUDA();
	void ComputeSolidForceCUDA();

	void computeNumBlocks(int numPnts, int maxThreads, int &numBlocks, int &numThreads);


	}
#endif