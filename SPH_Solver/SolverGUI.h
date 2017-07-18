#pragma once
#include <time.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>


#include "fluid_system_host.cuh"	
#include "fluid_system.h"

#include "stokesian/stokesian.h"

#include "CatToolBox.h"

struct vertex {
	cfloat3 pos;
	cfloat4 color;
};

typedef std::vector<vertex> varray;

class SolverGUI{

public:
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



	void Initialize(int argc,char** argv);
	
	void CreateCube(void);
	void DestroyCube(void);
	void DrawCube(void);
	void ReSize(int width,int height);

	
	void DrawParticleSystem();
	void CreateParticleSystem();



	//rendering particles
	GLuint ProjectionMatrixUniformLocation,
		ViewMatrixUniformLocation,
		ModelMatrixUniformLocation,
		ParticleSizeUniformLocation;
	GLuint ShaderIds[4]; //program, vertex, geometry, fragment
	GLuint BufferIds[3];

	//rendering simple geometry
	GLuint SimpleGeoShader[3]; //0-program, 1-vertex, 2-fragment
	varray planeGrid;

	cmat4 ProjectionMatrix,
		ViewMatrix,
		ModelMatrix;

	float CubeRotation = 0;
	clock_t LastTime = 0;

	cCamera camera;
	//buffer bindings
	displayPack* dispBuffer;
	int* pnum;





	//simulation
	FluidSystem* psys;
	stokesianSolver* stokesian;

	int maxPointNum;
	int pointNum;
};
