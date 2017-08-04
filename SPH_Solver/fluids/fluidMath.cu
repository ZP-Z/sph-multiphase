#include "fluidMath.cuh"
#include <stdio.h>
#include "geometry.h"

__device__ void multiply_matrix3(float* a, float* b, float* c){
	float d[9];
	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			d[i*3+j] = a[i*3+0]*b[0*3+j]+a[i*3+1]*b[1*3+j]+a[i*3+2]*b[2*3+j];
	for(int k=0; k<9; k++)
		c[k] = d[k];
}
__device__ __host__ void mat_x_vec(float* mat, float* vec, float* res){
	for(int i=0; i<3; i++){
		res[i] = 0;
		for(int j=0; j<3; j++){
			res[i]+=mat[i*3+j]*vec[j];
		}
	}
}

__device__ void transpose3(float* a,float* b){
	float c[9];
	c[0]=a[0]; c[1]=a[3]; c[2]=a[6];
	c[3]=a[1]; c[4]=a[4]; c[5]=a[7];
	c[6]=a[2]; c[7]=a[5]; c[8]=a[8];
	for(int k=0; k<9; k++)
		b[k]=c[k];
}
__device__ float det(float* a){
	float det = a[0]*a[4]*a[8] + a[1]*a[5]*a[6] + a[2]*a[3]*a[7];
	det -= (a[2]*a[4]*a[6] + a[1]*a[3]*a[8] + a[5]*a[7]*a[0]);
	return det;
}

