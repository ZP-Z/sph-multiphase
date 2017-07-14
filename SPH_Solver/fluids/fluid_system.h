

#pragma once

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common_header.h"
#include "../CatToolBox.h"

#define GRID_UCHAR			0xFF
#define GRID_UNDEF			4294967295

struct ParamCarrier {

	//environment
	float intstiff;
	float extstiff;
	float extdamp;
	float acclimit;
	float vlimit;
	cfloat3 gravity;
	cfloat3 volmin;
	cfloat3 volmax;
	cfloat3 softminx;
	cfloat3 softmaxx;

	//case
	int num;
	float dt;
	float simscale;

	//fluid
	float viscosity;
	float restdensity;
	float mass;
	float radius;
	float smoothradius;
	float scalep;
	float massArr[3];
	float densArr[3];
	float viscArr[3];


	//boundary particle
	float bmass;
	float bRestdensity;


	//grid parameter
	float cellsize;
	int searchnum;
	int neighbornum;
	cfloat3 gridmin,gridmax;
	cfloat3 gridsize;
	cfloat3 gridIdfac;
	cint3 gridres;
	int gridtotal;

	//tensile instability
	float w_deltax;
	float fijn;
	float a_phi;
	float a_psi;
	float k_c;

	float SFlipControl;
	float LFlipControl;
	int   mpmSplitNum;
	float collisionStiff;
	float surfaceTensionK;
	float surfaceTensionFluidC;
	float fBlending;
};

#define ROUND_EMIT 0
#define SQUARE_EMIT 1

class ParticleEmitter{
public:
	cfloat3 dir,
		position;
	float radius,
		spacing,
		vel;
	int timer,
		interval;
	void EmitImmediate(std::vector<cfloat3>& container);
};


struct displayPack {
	cfloat3 pos;
	cfloat4 color;
	int type;
};
struct calculationPack {
	cfloat3 vel;
	cfloat3 veleval;
	float pressure;
	float dens;
	float restdens;
	float restmass;
	float visc;
};
//remain unsorted
struct IntermediatePack {

};

struct bufList {
	//Particle properties
	cfloat3*			mpos;
	cfloat3*			mvel;
	cfloat3*			mveleval;
	cfloat3*			mforce;
	float*			mpress;
	float*			last_mpress;
	float*			mdensity;
	uint*			mgcell;
	uint*			mgndx;
	uint*			mclr;
	cfloat4*			mColor;
	int*			misbound;
	cfloat3*			accel; //leapfrog integration
							   //End particle properties

	displayPack* displayBuffer;
	calculationPack* calcBuffer;
	IntermediatePack* intmBuffer;

	//multi fluid
	float*			mf_alpha;				// MAX_FLUIDNUM * 4 bytes for each particle
	float*			mf_alpha_pre;			// MAX_FLUIDNUM * 4 bytes for each particle
	float*			mf_pressure_modify;	//  4 bytes for each particle
	cfloat3*			mf_vel_phrel;			// MAX_FLUIDNUM * 12 bytes for each particle
	float*			mf_restmass;
	float*			mf_restdensity;
	float*			mf_visc;
	cfloat3*			mf_velxcor;
	cfloat3*			mf_alphagrad;			// MAX_FLUIDNUM * 12 bytes for each particle

	int*            MFtype;
	int*			MFid;					//particle id
	int*			MFidTable;				//id table
	float*			MFtensor;				//sigma_s
	float*			MFtemptensor;			//sigma
	float*			MFRtensor;				//Artificial Tensor
	float*			MFvelgrad;				//velocity grad
	float*			MFpepsilon;

	int*			mf_multiFlag;			// basically a variously used buffer
											//End multi fluid

	uint*			mcluster;

	//For sorting
	char*			msortbuf;
	uint*			mgrid;
	int*			mgridcnt;
	int*			mgridoff;
	int*			mgridactive;
	//new sort
	uint*			midsort;
	//End sorting

	//Mpm Grid
	//Split Grid For Phases!
	float*          mpmMass;
	cfloat3*         mpmVel;
	cfloat3*         mpmForce;
	float*          mpmTensor;

	float*          mpmAlpha;

	uint*           mpmGid;    //mpmSize
	uint*           mpmIdSort; //mpmSize
	int*            mpmGridVList; //mpmSize
	int*            mpmGridCnt;//gridSize
	int*            mpmGridOff; //gridSize

	cfloat3*         mpmPos;

};// End particle&grid buffers

class FluidSystem {
public:
	FluidSystem ();
		
	


	// Setup
	void Setup ( bool bStart );
		
	void ResetParameters ();
	void ParseXML();
	void SetupSpacing ();
	void SetupGridAllocate ();
	void BeforeFirstRun();
	void SetupDevice();
	
	void SetupSimpleSphCase();


	//Manipulating particles
	void AllocateParticles(int cnt);
	int AddParticle();
	void InsertParticles(std::vector<cfloat3> plist,
		int cat);

	void saveParticle(std::string name);
	int loadParticle(std::string name);

	//Code Generate Entities
	void SetupMfAddVolume( cfloat3 min, cfloat3 max, float spacing, cfloat3 offs, int cat);// cat: category label
	//void SetupMfAddDeformVolume( cfloat3 min, cfloat3 max, float spacing, cfloat3 offs, int type);
	void SetupAddSolid( cfloat3 min, cfloat3 max, float spacing, int cat);
    void SetupAddBubble( cfloat3 min, cfloat3 max, float spacing, int constitution);

	//Code Generate Boundaries
    void LoadBoundary(std::string boundfile);
	void CalcBoundarySpacing();
    void SetupAddWall(cfloat3 min, cfloat3 max);
	void SetupAddOpenBox(cfloat3 min, cfloat3 max, float thickness);
    void SetupAddSphere(cfloat3 center, float radius);
    void SetupAddShell(cfloat3 center, float innner, float outter, float cutlevel);
    std::vector<cfloat3> boundaryPArray;

	//Code Emit Particles
	void EmitParticles(int cat);





	//Memery Transfer
	void TransferToCUDA();
	void TransferFromCUDA();
	void TransferFromCUDAForLoad();
	void EmitUpdateCUDA(int startnum, int endnum);
	void TransferToCUDA(bufList& fbuf);
	void TransferFromCUDA(bufList& fbuf);






	// Simulation
	void Run (int w, int h);
	//void RunSimulateMultiCUDAFull();
    //void RunMpmSolution();
	void RunSimpleSPH();
	//void RunSolid();

	void LiquidToBubble();

	void CaptureVideo (int width, int height);
	void record ( int param, std::string, cTime& start );
	void outputFile();
	
	void ClearTimer();
	void CheckTimer(const char* msg);



	void Exit();





	//VARIABLES

	// Particle Buffers
	//int						mNumPoints;
	//int						mMaxPoints;

	cfloat3*				mPos;
	DWORD*					mClr;
	cfloat4*				mColor;

	int*					mIsBound;
	cfloat3*				mVel;
	cfloat3*				mVelEval;
	float*					mPressure;
	float*					mDensity;
	cfloat3*				mForce;

	float*					m_alpha;  //size: mMaxPoints * MAX_FLUIDNUM
	float*					m_alpha_pre; //size: mMaxPoints * MAX_FLUIDNUM
	float*					m_pressure_modify; //size: mMaxPoints * MAX_FLUIDNUM
	cfloat3*				m_vel_phrel; //size: mMaxPoints * MAX_FLUIDNUM

	float*					m_restMass;
	float*					m_restDensity;
	float*					m_visc;
	cfloat3*				m_velxcor; //XSPH correction
	cfloat3*				m_alphagrad; //size: mMaxPoints * MAX_FLUIDNUM
		
	int*					MF_type;  //(project-u)
	float*					MF_tensor;//(project-u)
	int*					MF_id;
	int*					MF_idTable;
	float*					MF_pepsilon;
	
	displayPack* displayBuffer;
	calculationPack* calculationBuffer;
	IntermediatePack* intmBuffer;

	float boundarySpacing;
	float boundarySpacingFactor;

    int runMode,
		displayMode,
		sceneID,
		frameNo,
		outputNo;
	bool bPause,
		bSnapshot,
		bOutput;
	
	int pointNum;

	cfloat3 fvmin[3], fvmax[3];
	int fluidnum;
	int maxPointNum;
	float massratio[3];
	float densratio[3];
	float viscratio[3];
	float pSpacing, pRealDist;
	
	// Timer
	float lastTime;
	std::vector<const char*> timerStrs;
	std::vector<float> timerVals;

	ParamCarrier hostCarrier;
};