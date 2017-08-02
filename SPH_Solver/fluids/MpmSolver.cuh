#pragma once

#include "common_header.h"


__global__ void initMpm(bufList buf, int mpmSize);

__global__ void MpmColorTest(bufList buf, int mpmSize);

__global__ void MpmGetMomentum(bufList buf, int mpmSize);