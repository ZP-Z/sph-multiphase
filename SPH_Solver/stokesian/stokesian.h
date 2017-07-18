#pragma once

#include "Mobility.cuh"
#include "geometry.h"
#include <iostream>

//using namespace Eigen;
using namespace std;

__inline__ float randf(float max) {
	return (float)rand()/RAND_MAX * max;
}

class tensor3 {

public:
	tensor3() {data=NULL;}
	void setSize(int x, int y, int z){
		xl = x; yl = y; zl = z;
		if(data!=NULL)
			delete data;
		data = new float[xl*yl*zl];
		for (int i=0; i<xl*yl*zl; i++)
			data[i] = 0;
	}
	~tensor3() {
		if (data != NULL){
			delete data;
			data = NULL;
		}
	}

	float& operator() (int rx, int ry, int rz) {
		return data[rz*xl*yl + ry*xl + rx];
	}
	void printTensor() {
		for (int i=0; i<zl; i++) {
			printf("{");
			for (int j=0; j<yl; j++) {
				for (int k=0; k<xl; k++) {
					printf("%f ", (*this)(i, j, k));
				}
				if (j==yl-1)
					printf("}\n");
				else
					printf("\n");
			}
		}
	}
private:
	int xl, yl, zl;
	float* data;
};



class stokesianSolver {

public:

	//interpolation for pairwise resistence matrix
	tensor3 signature;
	//Matrix3d kroneckerDelta;

	int pointNum;
	int ndim; // 6n
	int ng;	  // 3n
	int nstep;
	float dt;
	int t;
	cfloat3 gravity;
    float alpha;
    float radius;

	float dr;

	displayPack* displayBuffer;
	cfloat3* direction;
	cfloat3* force; //external force
	cfloat3* torque;
	cfloat3* unew;
	cfloat3* omega; //angle speed
	//float* psize;
	float* lifetime;


    //cuda device buffer
	//cfloat3 * cuAMatrix;
    stokesianBufList deviceBuffer;
	

	//sorting grid
	//sortingGrid sGrid;
	/*float3 gxmin;
	float3 gxmax;
	float gcelldx;*/

	int frameNo;
	float simTimer;
	//paramReader pmReader;
	int particleCount;
	int objectCount;
	int maxAllowNum;
	
	int firstHalf;

	stokesianSolver();

	void SetupSolver();
	void LoadParam();
	//void readData();

	int pnum() {
		return pointNum;
	}
	int maxAllow(){
		return maxAllowNum;
	}

	//velocity superposition method
    //void SolveMobility();
	//void SolveMobilityGPU();
	//void SolveMobilityGPU_cutoff();
	
	void SetForce();
	void AllocateBuffers();

	
	void step();
	void SolveMobilityGPU_walkthrough();
	
	//run
	void RunMobilityGPU();

	void AddShape();
	void AddSphere(cfloat3 centre, float radius, float spacing,int id, float psize);
	//void AddBox(float3 xmin,float3 len,float spacing,int id);
	
	int BirthAndDeath();
	void ReleaseBubble();
	void InjectWater();


	~stokesianSolver() {
		/*if(pos)
			delete pos;*/
	}
};