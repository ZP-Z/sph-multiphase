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




class shaderObject {
public:
	GLuint programid;
	GLuint shaderid[3];//vertex, geometry, fragment
	shaderObject() {
		programid = glCreateProgram();
		shaderid[0] = 0;
		shaderid[1] = 0;
		shaderid[2] = 0;
	}
	
	void Release() {
		printf("A shader program is being destroyed.\n");
		for (int i=0; i<3; i++) {
			if (shaderid[i]>0) {
				glDetachShader(programid, shaderid[i]);
				glDeleteShader(shaderid[i]);
			}
		}
		glDeleteProgram(programid);
	}

	void loadShader(const char* filename, GLuint sdtype) {
		switch (sdtype) {
		case GL_VERTEX_SHADER:
			shaderid[0] = LoadShader(filename,sdtype);
			break;
		case GL_GEOMETRY_SHADER:
			shaderid[1] = LoadShader(filename,sdtype);
			break;
		case GL_FRAGMENT_SHADER:
			shaderid[2] = LoadShader(filename, sdtype);
			break;
		}
	}
	void LinkProgram() {
		for (int i=0; i<3; i++) {
			if (shaderid[i]>0) {
				glAttachShader(programid, shaderid[i]);
			}
		}
		glLinkProgram(programid);
		ExitOnGLError("ERROR: Could not link the shader program");
	}
};


#define PARTICLE_RENDERER 0
#define CUBE_RENDERER 1

class SolverGUI{

public:
	void SetupSolver();
	void GetBoundingBox();

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
	void InitializeGL(int argc,char** argv);
	void LoadParam();

	void ReleaseGLBuffers(void);
	
	void ReSize(int width,int height);

	void CreateShaders();
	void CreateCubicShaders();
	

	void DrawParticleSystem();
	void DrawParticleSystem_Cube();

	void AllocateParticleVBO();
	void AllocateCubeVBO();

	void AllocateBoxVBO();
	void CreateBoxShaders();
	void DrawBox();
	void TakeSnapshot();

	//rendering particles
	GLuint ProjectionMatrixUniformLocation,
		ViewMatrixUniformLocation,
		ModelMatrixUniformLocation,
		ParticleSizeUniformLocation;

	GLuint boxUniformLocation[3];//model-matrix, view-matrix, projection-matrix
		
	
	GLuint BufferIds[3];   //vao, vbo-(pos,color), vbo-(rotation)
	GLuint GeoBufferIds[3];//vao, vbo-(pos,color), indexbuffer

	//textures
	GLuint textures[10];


	//rendering simple geometry
	varray planeGrid;
	vertex boundingBox[8];


	cmat4 ProjectionMatrix,
		ViewMatrix,
		ModelMatrix;

	float CubeRotation = 0;
	clock_t LastTime = 0;

	cCamera camera;
	int shaderId;

	//buffer bindings
	displayPack* dispBuffer;
	int* pnum;
	cmat4* rotationBuffer;


	//shaders
	vector<shaderObject> shaders;


	//simulation
	FluidSystem* psys;
	stokesianSolver* stokesian;

	int maxPointNum;
	int pointNum;

	int rendermode;
	int frameNo;
	bool bSnapshot;
	
};
