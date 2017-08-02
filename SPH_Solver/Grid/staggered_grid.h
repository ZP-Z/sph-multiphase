#pragma once

#include "geometry.h"
#include <vector>

//base grid for flip and mpm
class Grid {
public:
	cint3 res;	//cell resolution
	cfloat3 xmin,xmax;
	float cellsize;

	int getId(cint3 idvec) {
		return idvec.y * (res.x*res.z) + idvec.z * res.x + idvec.x;
	}

	cint3 getIdvec(cfloat3 pos) {
		pos = (pos - xmin)/cellsize;
		return cint3(pos.x, pos.y, pos.z);
	}
	
};

class Staggered_Grid : public Grid {
public:
	cint3 ures;
	cint3 vres;
	cint3 wres;
	float *u, *v, *w;

	void SetRes(int x, int y, int z) {
		res.Set(x,y,z);
		ures.Set(x+1,y,z);
		vres.Set(x,y+1,z);
		wres.Set(x,y,z+1);

		u = new float[ ures.x*ures.y*ures.z ];
		v = new float[ vres.x*vres.y*vres.z ];
		w = new float[ wres.x*wres.y*wres.z ];
	}
	void Release() {
		delete u;
		delete v;
		delete w;
	}
};