

#pragma once

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common_header.h"
#include "../CatToolBox.h"

#define MAX_PARAM			60
#define GRID_UCHAR			0xFF
#define GRID_UNDEF			4294967295

//From YanXiao
#define OUTPUT_INT			1
#define START_OUTPUT		2
#define SHOW_BOUND			3
#define SAVE_STAT			4
#define CHANGE_DEN			6
#define SHOW_PROPERTY		7
#define SHOW_TYPE			8



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