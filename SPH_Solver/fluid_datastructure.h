#pragma once

#include "geometry.h"

struct ParamCarrier {

	//environment
	float intstiff;
	float extstiff;
	float extdamp;
	float acclimit;
	float vlimit;
	cfloat3 gravity;
	cfloat3 volmin;
	cfloat3 volmax;
	cfloat3 softminx;
	cfloat3 softmaxx;

	//case
	int num;
	float dt;
	float simscale;

	//fluid
	float viscosity;
	float restdensity;
	float mass;
	float radius;
	float smoothradius;
	float scalep;
	float massArr[3];
	float densArr[3];
	float viscArr[3];


	//boundary particle
	float bmass;
	float bRestdensity;


	//grid parameter
	float cellsize;
	int searchnum;
	int neighbornum;
	cfloat3 gridmin, gridmax;
	cfloat3 gridsize;
	cfloat3 gridIdfac;
	cint3 gridres;
	int gridtotal;

	//tensile instability
	float w_deltax;
	float fijn;
	float a_phi;
	float a_psi;
	float k_c;

	float SFlipControl;
	float LFlipControl;
	int   mpmSplitNum;
	float collisionStiff;
	float surfaceTensionK;
	float surfaceTensionFluidC;
	float fBlending;
};
