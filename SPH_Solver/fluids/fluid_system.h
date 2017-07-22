

#pragma once

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common_header.h"
#include "CatToolBox.h"

#define GRID_UCHAR			0xFF
#define GRID_UNDEF			4294967295

#define FLUID 0
#define BOUNDARY 1

#define ROUND_EMIT 0
#define SQUARE_EMIT 1

class ParticleEmitter {
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
	void Setup ();
		
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
	void InsertParticles(
		ParticleEmitter& emitter,
		std::vector<cfloat3> plist,
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
	//void TransferToCUDA();
	//void TransferFromCUDA();
	//void TransferFromCUDAForLoad();
	void EmitUpdateCUDA(int startnum, int endnum, bufList& fbuf);
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