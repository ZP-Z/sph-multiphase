#include <assert.h>
#include <stdio.h>
#include <conio.h>

//#include "gl_helper.h"
//#include <gl/glut.h>

//#include "camera3d.h"
//#include "common_defs.h"
//#include "mtime.h"
#include "fluid_system.h"
#include "fluid_system_host.cuh"


void FluidSystem::EmitUpdateCUDA (int startnum, int endnum)
{ 
	//int numpoints = endnum - startnum;
	//CopyEmitToCUDA ( (float*) mPos, (float*) mVel, (float*) mVelEval, (float*) mForce, mPressure, mDensity, NULL, NULL, (char*) mClr, startnum, numpoints , mIsBound); 

	//CopyEmitMfToCUDA ( m_alpha, m_alpha_pre, m_pressure_modify, (float*) m_vel_phrel, m_restMass, m_restDensity, m_visc, (float*)m_velxcor, (float*)m_alphagrad, startnum, numpoints);
	//CopyEmitToCUDA_Uproject((int*) MF_type, MF_tensor, MF_id, startnum, numpoints);

	//UpdatePNumCUDA(endnum);
	//cudaThreadSynchronize ();
}
