#pragma once

#include "common_header.h"

__global__ void ComputeSolidTensor(bufList fbuf, int pnum);
__global__ void ComputeSolidForce(bufList fbuf, int pnum);



//with reference shape
__device__ void contributeDeformGrad(cmat3& res, int i, bufList buf, int cell);
__global__ void ComputeSolidTensor_X(bufList buf, int pnum);
__device__ void contributeSolidForce_X(bufList buf, int pnum);
__global__ void ComputeSolidForce_X(bufList buf, int pnum);

__global__ void ComputeInvA(bufList buf,int pnum);