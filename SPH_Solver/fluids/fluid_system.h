

#pragma once

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../common_header.h"
#include "CatToolBox.h"
#include "fluid_system_kern.cuh"

#define GRID_UCHAR			0xFF
#define GRID_UNDEF			4294967295



#define ROUND_EMIT 0
#define SQUARE_EMIT 1


#define SPH 0
#define IISPH 1
#define SolidSPH 2
#define MPM 3

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
	void ParseXML(int caseid);
	void SetupSpacing ();
	void SetupGrid ();
	void SetupDevice();
	

	void SetupSimpleSphCase();
	void SetupSolidCase();


	void AllocateParticles(int cnt);
	int AddParticle();
	void InsertParticles(
		ParticleEmitter& emitter,
		std::vector<cfloat3> plist,
		int cat);

	void saveParticle(std::string name);
	int loadParticle(std::string name);




	//Code Generate Entities
	void AddFluidVolume( cfloat3 min, cfloat3 max, float spacing, cfloat3 offs, int cat);// cat: category label
	//void SetupMfAddDeformVolume( cfloat3 min, cfloat3 max, float spacing, cfloat3 offs, int type);
	void AddDeformableVolume( cfloat3 min, cfloat3 max, float spacing, int cat);

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
	void EmitUpdateCUDA(int startnum, int endnum, bufList& fbuf);
	void TransferToCUDA(bufList& fbuf);
	void TransferFromCUDA(bufList& fbuf);

	// Simulation
	void BeforeFirstStep();
	void Run();
	//void RunSimulateMultiCUDAFull();
    //void RunMpmSolution();
	void RunSimpleSPH();
	void RunSolid();
	void RunIISPH();
	void RunMPM();

	void LiquidToBubble();

	void outputFile();
	
	void ClearTimer();
	void CheckTimer(const char* msg);

	void Exit();



	//particle buffer
	displayPack* displayBuffer;
	calculationPack* calculationBuffer;

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

	


	//mpm 
	cfloat3 mpmxmin, mpmxmax;
	cint3 mpmres;
	float mpmcellsize;
	cfloat3* nodevel;
	
	void SetupMPMGrid();
	void ReleaseMPMGrid();
	void CopyMPMFromDevice();
};