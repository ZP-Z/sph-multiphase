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

void FluidSystem::TransferToCUDA ()
{ 
	CopyToCUDA ( (float*) mPos, (float*) mVel, (float*) mVelEval, (float*) mForce, mPressure, mDensity, NULL, NULL, (char*) mClr); 

	CopyMfToCUDA ( m_alpha, m_alpha_pre, m_pressure_modify, (float*) m_vel_phrel, m_restMass, m_restDensity, m_visc, (float*)m_velxcor, (float*)m_alphagrad);
	CopyToCUDA_Uproject((int*) MF_type, MF_tensor, MF_id);

	CopyBoundToCUDA(mIsBound);
}


void FluidSystem::TransferFromCUDA ()	
{
	CopyFromCUDA ( (float*) mPos, (float*) mVel, (float*) mVelEval, (float*) mForce, mPressure, mDensity, NULL, NULL, (char*) mClr, 1);
	CopyMfFromCUDA ( m_alpha, m_alpha_pre, m_pressure_modify, (float*) m_vel_phrel, m_restMass, m_restDensity, m_visc, (float*)m_velxcor, (float*)m_alphagrad, 1);
	CopyFromCUDA_Uproject(MF_type, MF_id, MF_pepsilon, MF_tensor, 1);

	CopyBoundFromCUDA(mIsBound);
}

void FluidSystem::LiquidToBubble() {
	CopyFromCUDA((float*)mPos, (float*)mVel, (float*)mVelEval, (float*)mForce, mPressure, mDensity, NULL, NULL, (char*)mClr, 1);
	CopyFromCUDA_Uproject(MF_type, MF_id, MF_pepsilon, MF_tensor, 1);

	for (int i=0; i< pointNum; i++) {
		if ( dot(mPos[i],mPos[i]) <100.0f ) {
			MF_type[i] = 2; //bubble
		}
	}

	//CopyToCUDA((float*)mPos, (float*)mVel, (float*)mVelEval, (float*)mForce, mPressure, mDensity, mClusterCell, mGridNext, (char*)mClr);
	CopyToCUDA_Uproject(MF_type, MF_tensor, MF_id);
}

void FluidSystem::TransferFromCUDAForLoad ()
{
	CopyFromCUDA ( (float*) mPos, (float*) mVel, (float*) mVelEval, (float*) mForce, mPressure, mDensity, NULL, NULL, (char*) mClr, 2);
	CopyMfFromCUDA ( m_alpha, m_alpha_pre, m_pressure_modify, (float*) m_vel_phrel, m_restMass, m_restDensity, m_visc, (float*)m_velxcor, (float*)m_alphagrad, 2);

	CopyBoundFromCUDA(mIsBound);
	CopyFromCUDA_Uproject(MF_type, MF_id, MF_pepsilon, MF_tensor, 2);
}

//void FluidSystem::prepareProju(){
//	param_proju[0]=coK;
//	param_proju[1]=coG;
//	param_proju[2]=phi;
//	param_proju[3]=coA;
//	param_proju[4]=coB;
//	param_proju[5]=coLambdaK;
//	param_proju[6]=cohesion;
//	param_proju[7]=boundaryVisc;
//	param_proju[8]=sleepvel;
//	param_proju[9]=m_Param[PSPACING] * m_Param[PSIMSCALE];
//
//	param_proju[10]=coN;
//	param_proju[11]=Yradius;
//	param_proju[12]=vfactor;
//	param_proju[13]=fpfactor;
//	param_proju[14]=spfactor;
//	param_proju[15]=fsA;
//	param_proju[16]=fsB;
//	param_proju[17]=bdamp;
//	param_proju[18]=coD;
//	param_proju[19]=coD0;
//
//	//Solid for Deformable, Non-prefix for Sand
//	param_proju[20]=solid_coG;
//	param_proju[21]=solid_coV;
//	param_proju[22]=solid_coK;
//	param_proju[23]=solid_coA; 
//	param_proju[24]=solid_coB;
//	param_proju[25]=solid_fsa;
//	param_proju[26]=solid_fsb;
//	param_proju[27]=solid_coN;
//	param_proju[28]=solid_phi;
//	param_proju[29]=solid_Yradius;
//	param_proju[30]=fluidVConstraint;
//	param_proju[31]=tohydro;
//    //32 is mpmSpacing
//    param_proju[33]=m_Vec[PVOLMIN].x;
//    param_proju[34]=m_Vec[PVOLMIN].y;
//    param_proju[35]=m_Vec[PVOLMIN].z;
//    param_proju[36]=m_Vec[PVOLMAX].x;
//    param_proju[37]=m_Vec[PVOLMAX].y;
//    param_proju[38]=m_Vec[PVOLMAX].z;
//}

void FluidSystem::EmitUpdateCUDA (int startnum, int endnum)
{ 
	int numpoints = endnum - startnum;
	CopyEmitToCUDA ( (float*) mPos, (float*) mVel, (float*) mVelEval, (float*) mForce, mPressure, mDensity, NULL, NULL, (char*) mClr, startnum, numpoints , mIsBound); 

	CopyEmitMfToCUDA ( m_alpha, m_alpha_pre, m_pressure_modify, (float*) m_vel_phrel, m_restMass, m_restDensity, m_visc, (float*)m_velxcor, (float*)m_alphagrad, startnum, numpoints);
	CopyEmitToCUDA_Uproject((int*) MF_type, MF_tensor, MF_id, startnum, numpoints);

	UpdatePNumCUDA(endnum);
	cudaThreadSynchronize ();
}
