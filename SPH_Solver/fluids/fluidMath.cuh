#pragma once

//#include <cutil.h>					// cutil32.lib
//#include <cutil_math.h>				// cutil32.lib
#include "host_defines.h"
#include "svd3_cuda.h"


__device__ void multiply_matrix3(float* a, float* b, float* c);
__device__ void transpose3(float* a,float* b);
__device__ float det(float* a);

__device__ __host__ void mat_x_vec(float* mat, float* vec, float* res);

