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

	#include <vector_types.h>	
	#include <driver_types.h>			// for cudaStream_t
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

	//void FluidSetupCUDA (  int num, int gsrch, int3 res, cfloat3 size, cfloat3 delta, cfloat3 gmin, cfloat3 gmax, int total, int chk);
	//void FluidParamCUDA ( float ss, float sr, float pr, float mass, float rest, cfloat3 bmin, cfloat3 bmax, float estiff, float istiff, float pbstiff, float visc, float damp, float fmin, float fmax, float ffreq, float gslope, float gx, float gy, float gz, float al, float vl );
	/*void FluidParamCUDA_projectu(float coK, float coG, float phi,float coA,
		float coB, float coLambdaK, float boundaryVisc, float sleepvel, float cohesion, 
		float initspacing, float coN, float Yradius, float v_factor, float f_factor, float s_factor, float fsA, float fsB, float bdamp,
		float coD,float coD0,
		float solid_coG, float solid_coV, float solid_coK, float solid_coA, float solid_coB, float solid_fsa, float solid_fsb, float solid_coN, float solid_phi, float solid_Yradius,float fluidVConstraint,float tohydro);
	*/
	void FluidParamCUDA(ParamCarrier& param);

	void FluidParamCUDAbuffer_projectu(float* buffer);
	

    void MpmAllocateBuffer();

	//multi fluid
	void FluidMfParamCUDA ( float *dens, float *visc, float diffusion, float catnum, float dt, cfloat3 cont, cfloat3 mb1, cfloat3 mb2, float relax,int example);
	void CopyMfToCUDA ( float* alpha, float* alpha_pre, float* pressure_modify, float* vel_phrel, float* restmass, float* restdensity, float* visc, float* velxcor, float* alphagrad);
	void CopyMfFromCUDA ( float* alpha, float* alpha_pre, float* pressure_modify, float* vel_phrel, float* restmass, float* restdensity, float* visc, float* velxcor, float* alphagrad, int mode);

	//emit
	void CopyEmitToCUDA ( float* pos, float* vel, float* veleval, float* force, float* pressure, float* density, uint* cluster, uint* gnext, char* clr, int startnum, int numcount ,int* mIsBound);
	void CopyEmitMfToCUDA ( float* alpha, float* alpha_pre, float* pressure_modify, float* vel_phrel, float* restmass, float* restdensity, float* visc, float* velxcor, float* alphagrad, int startnum, int numcount);
	void CopyEmitToCUDA_Uproject(int* mftype, float* tensorbuffer, int* bornid, int startnum, int numcount);
	void UpdatePNumCUDA( int newPnum);

	void CopyToCUDA ( float* pos, float* vel, float* veleval, float* force, float* pressure, float* density, uint* cluster, uint* gnext, char* clr );
	void CopyFromCUDA ( float* pos, float* vel, float* veleval, float* force, float* pressure, float* density, uint* cluster, uint* gnext, char* clr, int mode);

	void CopyBoundToCUDA(int* isbound);
	void CopyBoundFromCUDA(int* isbound);

	void CopyToCUDA_Uproject(int* mftype, float* tensorbuffer,int* bornid);
	void CopyFromCUDA_Uproject(int* mftype, int* idtable, float* pepsilon, float* stensor, int mode);

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
	void MfChangeDensityCUDA(const float scale);

	//Project-U changes computing force
	void ComputeForceCUDA_ProjectU();

    //Mpm Force
    void ComputeMpmForce();

	void ComputeSolidTensor();
	void ComputeSolidForce();
	}
#endif