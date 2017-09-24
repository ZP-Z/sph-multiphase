#pragma once

#include "fluid_system_kern.cuh"


__global__ void initMpm(bufList buf, int mpmSize);

__global__ void MpmColorTest(bufList buf, int mpmSize);

__global__ void MpmGetMomentum(bufList buf, int mpmSize);

__global__ void MpmParticleToGrid(bufList buf, int mpmSize);
__global__ void MpmParticleStress(bufList buf,int pnum);
__global__ void MpmParticleStressFiniteStrain(bufList buf, int pnum);
__global__ void MpmNodeUpdate(bufList buf,int mpmSize);
__global__ void MpmParticleUpdate(bufList buf,int pnum);

__global__ void MpmParticleDp(bufList buf, int pnum);
__global__ void MpmParticleToGrid_APIC(bufList buf,int mpmSize);
__global__ void MpmParticleUpdate_APIC(bufList buf,int pnum);

