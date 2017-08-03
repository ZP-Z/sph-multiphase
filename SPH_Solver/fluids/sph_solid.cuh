#pragma once

#include "common_header.h"

__global__ void ComputeSolidTensor(bufList fbuf, int pnum);
__global__ void ComputeSolidForce(bufList fbuf, int pnum);
