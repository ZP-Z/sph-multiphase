#pragma once

#include <math.h>
#include <stdio.h>

#include "GL/glew.h"
#include "GL/freeglut.h"
#include "SOIL.h"
#include "CatTimer.h"

#ifdef BUILD_CUDA
#include <host_defines.h>
#define HDFUNC __host__ __device__
#endif

HDFUNC inline float cmin(float a, float b) {
	return a < b ? a : b;
}


inline float cmax(float a, float b) {
	return a > b ? a : b;
}

struct cfloat2 {
	float x,y;
	cfloat2(){}
	cfloat2(float _x,float _y):x(_x),y(_y){}
	void Set(float _x,float _y){ x=_x; y=_y; }
	cfloat2 operator*(float s){
		return cfloat2(x*s, y*s);
	}
	cfloat2 operator/(float s){
		return cfloat2(x/s,y/s);
	}
	cfloat2 operator-(cfloat2& b){
		return cfloat2(x-b.x,y-b.y);
	}
	cfloat2 operator+(cfloat2& b){
		return cfloat2(x+b.x, y+b.y);
	}
};

inline float dot(cfloat2& a, cfloat2& b) {
	return a.x*b.x + a.y*b.y;
}


struct cfloat3 {
	float x, y, z;
	HDFUNC cfloat3() {}
	HDFUNC cfloat3(float _x, float _y, float _z) :x(_x), y(_y), z(_z) {}
	HDFUNC void Set(float _x, float _y, float _z) {
		x = _x;
		y = _y;
		z = _z;
	}
	HDFUNC cfloat3 operator- (cfloat3& b) {
		return cfloat3(x - b.x, y - b.y, z - b.z);
	}
	HDFUNC cfloat3 operator+(cfloat3& b) {
		return cfloat3(x + b.x, y + b.y, z + b.z);
	}
	HDFUNC cfloat3 operator*(cfloat3& b) {
		return cfloat3(x * b.x, y * b.y, z * b.z);
	}
	HDFUNC cfloat3 operator/(cfloat3& b) {
		return cfloat3(x / b.x, y / b.y, z / b.z);
	}
	HDFUNC cfloat3& operator/(float b){
		return cfloat3(x/b,y/b,z/b);
	}

	HDFUNC cfloat3 operator*(float s) {
		return cfloat3(x*s, y*s, z*s);
	}
	HDFUNC cfloat3 operator+(float a) {
		return cfloat3(x + a, y + a, z + a);
	}
	HDFUNC cfloat3 operator-(float a) {
		return cfloat3(x - a, y - a, z - a);
	}
	HDFUNC void operator+= (cfloat3& b) {
		x+=b.x; y+=b.y; z+=b.z;
	}
	HDFUNC void operator-= (cfloat3& b){
		x-=b.x; y-=b.y; z-=b.z;
	}
	HDFUNC void operator*= (float b) {
		x*=b; y*=b; z*=b;
	}
	HDFUNC void operator/= (float b) {
		x/=b; y/=b; z/=b;
	}
	HDFUNC void operator /= (cfloat3& b){
		x/=b.x; y/=b.y; z/=b.z;
	}

	HDFUNC float minx() {
		return cmin(x, cmin(y, z));
	}
};

HDFUNC inline cfloat3 minfilter(cfloat3 a, cfloat3 b) {
	return cfloat3(cmin(a.x, b.x), cmin(a.y, b.y), cmin(a.z, b.z));
}

inline cfloat3 maxfilter(cfloat3 a, cfloat3 b) {
	return cfloat3(cmax(a.x, b.x), cmax(a.y, b.y), cmax(a.z, b.z));
}

inline cfloat3 cross(cfloat3& a,cfloat3& b);
HDFUNC inline float dot(cfloat3& a,cfloat3& b){
	return a.x*b.x+a.y*b.y+a.z*b.z;
}

float angle(cfloat3& a,cfloat3& b);

struct cint3 {
	int x, y, z;
	HDFUNC cint3() {}
	HDFUNC cint3(int _x, int _y, int _z) :x(_x), y(_y), z(_z) {}
	HDFUNC void Set(int _x, int _y, int _z) {
		x = _x;
		y = _y;
		z = _z;
	}
	HDFUNC operator cfloat3(){
		return cfloat3(x,y,z);
	}
};



struct cmat3 {
	float data[9];
	cmat3() {}
	void Set(float mat[]) {
		for (int i=0; i<9; i++) data[i] = mat[i];
	}
	void Print() {
		for (int i=0; i<3; i++) {
			printf("%f %f %f\n", data[i*3], data[i*3+1], data[i*3+2]);
		}
	}
	inline float* operator[](int i){
		return &data[i*3];
	}
};

void mat3add(cmat3& a, cmat3& b, cmat3& c);
void mat3sub(cmat3& a, cmat3& b, cmat3& c);
void mat3prod(cmat3& a, cmat3& b, cmat3& c);

void mvprod(cmat3& m, cfloat3& v, cfloat3& c);

void RotateX(cfloat3& a, float b);
void RotateY(cfloat3& a, float b);
void RotateZ(cfloat3& a, float b);
void RotateXYZ(cfloat3& a, cfloat3& xyz);

class Grid {
protected:
	float cellsize;
	cint3 res;
	cfloat3 xmin, xmax;
	float* data;

};





/* -----------------------------------


	            OpenGL Utility


-------------------------------------*/

const float cPI = 3.1415926536;

struct cvertex{
	float position[4];
	float color[4];
};

struct cfloat4{
	float x,y,z,w;
	void Set(float _x,float _y,float _z,float _w){
		x=_x; y=_y; z=_z; w=_w;
	}
	cfloat4():
		x(0),y(0),z(0),w(0){}
	cfloat4(float _x, float _y, float _z, float _w):
		x(_x),y(_y),z(_z),w(_w){}
};

//watch out, column major
struct cmat4{
	float data[16];
	float* operator[](int i){
		return &data[i*4];
	}
	cmat4 operator*(cmat4& b){
		cmat4 tmp;
		for(int i=0;i<4;i++){
			for(int j=0;j<4;j++){
				tmp[i][j] = data[0+i*4]*b[0][j] + data[1+i*4]*b[1][j]
					+ data[2+i*4]*b[2][j] + data[3+i*4]*b[3][j];
			}
		}
		return tmp;
	}
};


extern const cmat4 IDENTITY_MAT;

float cotangent(float angle);
float deg2rad(float deg);
float rad2deg(float rad);

void RotateAboutX(cmat4& m, float ang);
void RotateAboutY(cmat4& m, float ang);
void RotateAboutZ(cmat4& m, float ang);
void ScaleMatrix(cmat4& m, cfloat3 x);
void TranslateMatrix(cmat4& m, cfloat3 x);

cmat4 CreateProjectionMatrix(float fovy, float aspect_ratio,
	float near_plane, float far_plane);

typedef 
	unsigned char
 uchar;

class keystack{
public:
	unsigned char buffer[32];
	unsigned char maxlen,lenth;
	keystack(){
		maxlen = 32;
		lenth = 0;
		memset(buffer, 0, sizeof(buffer));
	}
	int findkey(uchar key){
		int tmp = -1;
		for (int i=0; i<lenth; i ++) {
			if (buffer[i]==key) {
				return i;
			}
		}
		return tmp;
	}
	void push(uchar key){
		if(lenth==31){
			printf("key stack full\n");
			return;
		}
		int id = findkey(key);
		if (id>=0)
			printf("key exists.\n");
		else{
			buffer[lenth] = key;
			lenth += 1;
		}	
	}

	void pop(uchar key){
		int id = findkey(key);
		if(id<0)
			printf("key to pop unfound.\n");
		else{
			for(int i=id; i<lenth-1;i++){
				buffer[i] = buffer[i+1];
			}
			lenth -= 1;
		}
	}
};


class cCamera{
public:

	cmat4 projmat;
	cmat4 viewmat;

	//projection matrix
	float fovy, aspect_ratio, nearclip, farclip;

	cfloat3 dir,up,right,forceup;
	cfloat3 pos, target;
	
	//keyboard translate
	cfloat3 vel;
	float velf, velMax;
	keystack wasd;

	//drag rotation
	cfloat2 rotatexy;
	cfloat2 transxy;
	//drag zoom
	float zoomin;
	
	cCamera(){
		velf = 0.0f;
		rotatexy.Set(0,0);
		transxy.Set(0,0);
	}

	//ViewMatrix Utility
	cmat4 CameraRotateMat();
	void CameraTranslateMat(cmat4& rotation, cfloat3& pos);
	
	//ProjectionMatrix Utility
	void SetProjParam(float f,float a,float nc,float fc){
		fovy=f; aspect_ratio=a; nearclip=nc; farclip=fc;
	}
	void ProjectionMat(){
		projmat = CreateProjectionMatrix(fovy, aspect_ratio, nearclip, farclip);
	}


	void lookat(cfloat3& _pos, cfloat3& _target);

	//Moving Utility
	void MoveCamera(cfloat3 dx){
		pos = pos+dx;
		target = target+dx;
		lookat(pos, target);
	}

#define FRONT_VEL 0
#define RIGHT_VEL 1
#define UP_VEL 2

	void getVelDir(){
		if(wasd.lenth==0)
			return;

		switch (wasd.buffer[wasd.lenth-1]) {
		case 'w': case 'W':
			vel = dir * (-1);
			break;
		case 's': case 'S':
			vel = dir ;
			break;
		case 'd': case 'D':
			vel = right ;
			break;
		case 'a': case 'A':
			vel = right * (-1);
			break;
		case 'q': case 'Q':
			vel = up;
			break;
		case 'z': case 'Z':
			vel = up * (-1);
			break;
		}
	}

	void SetVelocty(unsigned char key){
		key = tolower(key);
		wasd.push(key);
		
	}

	void UnSetVelocity(unsigned char key){
		key = tolower(key);
		wasd.pop(key);
	}

	void AdvanceCamera(float dt);


};


