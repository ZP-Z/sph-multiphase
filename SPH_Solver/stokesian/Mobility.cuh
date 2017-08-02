
#pragma once
//#include <cutil_math.h>
#include <cuda_runtime.h>
#include "geometry.h"
#include "common_header.h"
//#include "stokesian\NeighborGrid.cuh"

struct stokesianBufList {
	cfloat3 * cuUnew;
	displayPack* dispBuffer;
	cfloat3 * cuOmega;
	cfloat3 * cuForce;

	//params
	float dt;
	float simscale;
};

//set epsilonIJK signature
void setSignature();


void getMobU(float* mat,float* f,float* u, int nsize);

//void getMobU_cutoff(float3* pos, int nsize, float3* f, float3* u, sortingGrid grid);


void getMobU_walkthrough(stokesianBufList buflist, int pnum);