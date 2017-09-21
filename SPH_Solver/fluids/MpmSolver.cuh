#pragma once

#include "common_header.h"


__global__ void initMpm(bufList buf, int mpmSize);

__global__ void MpmColorTest(bufList buf, int mpmSize);

__global__ void MpmGetMomentum(bufList buf, int mpmSize);

__global__ void MpmParticleToGrid(bufList buf, int mpmSize);
__global__ void MpmParticleStress(bufList buf,int pnum);
__global__ void MpmNodeUpdate(bufList buf,int mpmSize);
__global__ void MpmParticleUpdate(bufList buf,int pnum);

