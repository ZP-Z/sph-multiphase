
//Core of the stokesian solver
//#include <cutil_math.h>
#include "Mobility.cuh"
#include "stokesian.h"


//Allocate Some Buffer
stokesianSolver::stokesianSolver() {
}

void stokesianSolver::SetupSolver() {
	LoadParam();
	//SetupParameters("stokesian/configuration.txt");
	AllocateBuffers();
	//sGrid.SetupGridBuffer(gxmin,gxmax,gcelldx,pointNum);
	pointNum=0;
	frameNo = 0;
	//AddShape();
	
}

void stokesianSolver::AllocateBuffers() {
	//Note: F-T version stokesian dynamics
	
	//host buffer
	displayBuffer       = new displayPack[maxAllowNum];
	direction = new cfloat3[maxAllowNum];
	rotationMat = new cmat4[maxAllowNum];
	force     = new cfloat3[maxAllowNum];
	//torque    = new cfloat3[maxAllowNum];
	unew      = new cfloat3[maxAllowNum];
	omega     = new cfloat3[maxAllowNum];
	//psize     = new float[maxAllowNum];

	//device buffers
	cudaMalloc(&deviceBuffer.cuUnew,  sizeof(cfloat3)*maxAllowNum);
	cudaMalloc(&deviceBuffer.cuOmega, sizeof(cfloat3)*maxAllowNum);
	cudaMalloc(&deviceBuffer.dispBuffer,   sizeof(displayPack)*maxAllowNum);
	cudaMalloc(&deviceBuffer.cuForce, sizeof(cfloat3)*maxAllowNum);

	//initialize
	setSignature();
}





#ifdef BUBBLE
void stokesianSolver::ReleaseBubble(){
	//Sample And Give Birth
	cfloat3 startpoint;
	int startnum = realpnum;

	if (frameNo % 10 == 0)
		for (int i=0; i<firstHalf; i++) {
			startpoint = displayBuffer[i].pos;
			pos[realpnum] = startpoint + cfloat3(randf(2), randf(2), randf(2));
			psize[realpnum] = -10; //anti-gravity
			realpnum ++;
		}

	//Kill Particles Out of Date
	pointNum = realpnum;
}
#endif

void stokesianSolver::InjectWater(){

	if (pointNum == maxAllowNum) {
		return;
	}

	if (frameNo % 22 != 0)
		return;
	//printf("inject at turn %d\n", frameNo);
	float emitradius = 5;
	float spacing = 2;

	int res = emitradius/spacing;

	int startnum = pointNum;

	for (int i=-res; i<=res; i++) {
		for (int j=-res; j<=res; j++) {
			
			if (pointNum == maxAllowNum) {
				break;
			}

			int x = i*spacing;
			int z = j*spacing;
			cfloat3 rand3((float)rand()/RAND_MAX, (float)rand()/RAND_MAX,(float)rand()/RAND_MAX);
			rand3 = rand3*spacing/4 - spacing/8;

			displayBuffer[pointNum].pos.Set(x, 50, z);
			displayBuffer[pointNum].pos += rand3;

			displayBuffer[pointNum].color.Set(1, 0, 0, 1);
			displayBuffer[pointNum].type = 0;

			direction[pointNum].Set(0, 1, 0);
			rotationMat[pointNum] = IDENTITY_MAT;

			//psize[pointNum] = 1;
			pointNum++;
			if (pointNum == maxAllowNum) {
				printf("maximum number reached\n");
			}

		}
	}
	if (pointNum > startnum) {
		cudaMemcpy( deviceBuffer.dispBuffer+startnum,  displayBuffer+startnum, (pointNum-startnum)*sizeof(displayPack), cudaMemcpyHostToDevice);
	}

}

int stokesianSolver::BirthAndDeath(){
	//ReleaseBubble();
	InjectWater();
	return 0;
}


bool tmp = false;
//Run the simulation

void stokesianSolver::step() {
	SetForce();
	RunMobilityGPU();
	frameNo++;
}

void rotateF3withMat(cmat4& rotationmat, cfloat3& a) {
	cfloat3 res;
	res.x = rotationmat[0][0]*a.x + rotationmat[1][0]*a.y + rotationmat[2][0]*a.z;
	res.y = rotationmat[0][1]*a.x + rotationmat[1][1]*a.y + rotationmat[2][1]*a.z;
	res.z = rotationmat[0][2]*a.x + rotationmat[1][2]*a.y + rotationmat[2][2]*a.z;
	a = res;
}

void stokesianSolver::RunMobilityGPU() {
	
	simTimer += dt;
	
	SolveMobilityGPU_walkthrough();

	float rotation[9];
	cfloat3 tmpf3;
	
	//Advance
	for (int i=0; i<pointNum; i++) {
		
		//rotation
		tmpf3 = omega[i]*dt;
		
		cfloat3 init(0, 1, 0);
		//RotateXYZ(direction[i], tmpf3);

		RotateAboutX(rotationMat[i], tmpf3.x);
		RotateAboutY(rotationMat[i], tmpf3.y);
		RotateAboutZ(rotationMat[i], tmpf3.z);
		rotateF3withMat(rotationMat[i], init);
		
		//printf("%f %f %f %f %d\n",dot(direction[i],direction[i]),direction[i].x,direction[i].y,direction[i].z,i);
		displayBuffer[i].color.Set(fabs(init.x), fabs(init.y), fabs(init.z),1);
	}
	BirthAndDeath();
	
	//tmp = true;
	cudaMemcpy(deviceBuffer.dispBuffer, displayBuffer, pointNum*sizeof(displayPack), cudaMemcpyHostToDevice);
	
}

//Velocity Superposition Method

#define ROUND 0
#define SQUARE 1
#define SPHERE 2
#define BOX 3

#include <stdlib.h>

void stokesianSolver::AddShape() {
	
	int shape = 2;
	
	//round
	switch (shape) {
	case ROUND:{
		float rad = 20;
		float PI = 3.14159;
		for (int i=0; i<maxAllowNum; i++) {
			displayBuffer[i].pos.Set(rad * cos(PI/(maxAllowNum/2) * (i+1)),
				rad * sin(PI/(maxAllowNum/2) * (i+1)),
				1000);
		}
		break;
	}
	case SQUARE:{
		int roughsqrt = (int)sqrt((float)maxAllowNum)+1;
		float spacing = 8;
		float xmin = 0 - roughsqrt*spacing/2;
		float ymin = 0 - roughsqrt*spacing/2;
		
		for (int i=0; i<maxAllowNum; i++) {
			int x = i%roughsqrt;
			int y = i/roughsqrt;
			displayBuffer[i].pos.Set(xmin + x*spacing, ymin + y*spacing, 1500);
		}
		break;
	}
	case SPHERE: {
		int roughtroot = (int)pow((float)maxAllowNum,1/3.0);
		float spacing = 16;
		float sphRadius = roughtroot*spacing/2;
		float x,y,z;
		int count = 0;
		for (int i=0; i<roughtroot; i++) {
			for (int j=0; j<roughtroot; j++) {
				for (int k=0; k<roughtroot; k++) {
					z = i*spacing - sphRadius + ((float)rand()/RAND_MAX-0.5)/10*spacing;
					y = j*spacing - sphRadius + ((float)rand()/RAND_MAX-0.5)/10*spacing;
					x = k*spacing - sphRadius + ((float)rand()/RAND_MAX-0.5)/10*spacing;
					if(x*x + y*y + z*z >= sphRadius*sphRadius)
						continue;
					displayBuffer[pointNum].pos.Set(x, y, z);
					displayBuffer[pointNum].color.Set(1,0,0,1);
					pointNum++;
				}
			}
		}
		printf("sphere radius %f pnum %d\n",sphRadius, pointNum);
		break;
	}
	case BOX: {
		float xlen = 1;
		float ylen = 1;
		float zlen = 1;
		float spacing = 32;
		float ratio = pow(maxAllowNum/xlen/ylen/zlen,1/3.0);
		int xn = xlen*ratio;
		int yn = ylen*ratio;
		int zn = zlen*ratio;
		maxAllowNum = xn*yn*zn;
		int count = 0;
		for (int i=0; i<zn; i++) {
			for(int j=0;j<yn;j++){
				for(int k=0;k<xn;k++){
					displayBuffer[pointNum].pos.Set((k-xn/2) * spacing, (j-yn/2) * spacing,  (i-zn/2) * spacing+100);
					pointNum++;
				}
			}
		}
		printf("cube pnum %d\n",pointNum);
		break;
	}
		
	}
	
	//square
	
	
}

#include "tinyxml2.h"
#include "catXMLhelper.h"
using namespace tinyxml2;
extern XMLElement* pele;

void stokesianSolver::LoadParam() {
	tinyxml2::XMLDocument doc;
	int result = doc.LoadFile("stokesian/configuration.xml");
	printf("return by open xml %d\n", result);

	XMLElement* paramElement = doc.FirstChildElement("param");
	if (!paramElement)
	{
		printf("missing param xml node");
		exit(-1);
	}
	pele = paramElement;
	maxAllowNum = XMLGetFloat("pnum");
	dt = XMLGetFloat("dt");
	gravity=XMLGetFloat3("gravity");
	deviceBuffer.simscale = XMLGetFloat("simscale");
	deviceBuffer.dt = dt;

	pele = NULL;
}




void stokesianSolver::SetForce() {
	
	for (int i=0; i<pointNum; i++) {
		force[i]=gravity;
		//force[i].Set(0,0,0);

		//boundary
		if (displayBuffer[i].pos.y < 0)
			force[i].y += 0.5*(-(displayBuffer[i].pos.y));
		if (displayBuffer[i].pos.x < -30)
			force[i].x += 0.5*(-(displayBuffer[i].pos.x+30));
		if (displayBuffer[i].pos.x > 30)
			force[i].x += 0.5*(-(displayBuffer[i].pos.x-30));
		if (displayBuffer[i].pos.z < -30)
			force[i].z += 0.5*(-(displayBuffer[i].pos.z+30));
		if (displayBuffer[i].pos.z > 30)
			force[i].z += 0.5*(-(displayBuffer[i].pos.z-30));
	}

	cudaMemcpy(deviceBuffer.cuForce, force, sizeof(cfloat3)*pointNum, cudaMemcpyHostToDevice);
}

void stokesianSolver::SolveMobilityGPU_walkthrough() {
	getMobU_walkthrough(deviceBuffer, pointNum);
	
	cudaMemcpy(displayBuffer, deviceBuffer.dispBuffer, pointNum*sizeof(displayPack), cudaMemcpyDeviceToHost);
	cudaMemcpy(unew, deviceBuffer.cuUnew, pointNum*sizeof(cfloat3), cudaMemcpyDeviceToHost);
	cudaMemcpy(omega, deviceBuffer.cuOmega, pointNum*sizeof(cfloat3),cudaMemcpyDeviceToHost);
}

