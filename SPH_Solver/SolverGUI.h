#pragma once
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>

#include "fluid_system_host.cuh"	

#include "fluid_system.h"
#include "CatToolBox.h"


class SolverGUI{

public:

	//FluidSystem* psys;

	void SetupSolver();

	void Run();

	void Exit();

	//keyboard function
	void keyUp(unsigned char key);
	void keyDown(unsigned char key);

	//mouse function
	void drag(int dx,int dy);
	void MouseClick(int x,int y,int state, int button);

	//GL functions
	void render();

	//rendering utilities
	GLuint ProjectionMatrixUniformLocation,
		ViewMatrixUniformLocation,
		ModelMatrixUniformLocation,
		ParticleSizeUniformLocation,
		BufferIds[3]={0},
		ShaderIds[4]={0};
	cmat4 ProjectionMatrix,
		ViewMatrix,
		ModelMatrix;
	float CubeRotation = 0;
	clock_t LastTime = 0;

	cCamera camera;

	void Initialize(int argc,char** argv);
	
	void CreateCube(void);
	void DestroyCube(void);
	void DrawCube(void);
	void ReSize(int width,int height);

	int maxPointNum;
	int pointNum;
	void DrawParticleSystem();
	void CreateParticleSystem();

	//simulation
	FluidSystem* psys;
};


void cudaInit(int argc, char** argv);

