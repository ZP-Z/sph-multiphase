#include "geometry.h"



void  mat3add(cmat3& a, cmat3& b, cmat3& c) {
	for (int i=0; i<9; i++)
		c.data[i] = a.data[i]+b.data[i];
}

void  mat3sub(cmat3& a, cmat3& b, cmat3& c) {
	for (int i=0; i<9; i++)
		c.data[i] = a.data[i]-b.data[i];
}

void  mat3prod(cmat3& a, cmat3& b, cmat3& c) {
	for (int i=0; i<3; i++) {
		for (int j=0; j<3; j++) {
			c.data[i*3+j] = a.data[i*3]*b.data[j]+ a.data[i*3+1]*b.data[3+j] + a.data[i*3+2]*b.data[6+j];
		}
	}
}

void  mvprod(cmat3& m, cfloat3& v, cfloat3& c) {
	cfloat3 tmp;
	tmp.x = m.data[0]*v.x + m.data[1]*v.y + m.data[2]*v.z;
	tmp.y = m.data[3]*v.x + m.data[4]*v.y + m.data[5]*v.z;
	tmp.z = m.data[6]*v.x + m.data[7]*v.y + m.data[8]*v.z;
	c=tmp;
}


void  RotateX(cfloat3& a, float b) {
	float tmp[9]={1,0,0,    0,cos(b),-sin(b),    0,sin(b),cos(b)};
	cmat3 tmpmat;
	tmpmat.Set(tmp);
	mvprod(tmpmat, a, a);
}

void  RotateY(cfloat3& a, float b) {
	float tmp[9] ={cos(b),0,sin(b),    0,1,0,    -sin(b),0,cos(b)};
	cmat3 tmpmat;
	tmpmat.Set(tmp);
	mvprod(tmpmat, a, a);
}

void  RotateZ(cfloat3& a, float b) {
	float tmp[9] ={cos(b),-sin(b),0,    sin(b),cos(b),0,   0,0,1};
	cmat3 tmpmat;
	tmpmat.Set(tmp);
	mvprod(tmpmat, a, a);
}

void  RotateXYZ(cfloat3& a, cfloat3& xyz) {
	cfloat3 tmp = a; //backup
	RotateX(tmp, xyz.x);
	RotateY(tmp, xyz.y);
	RotateZ(tmp, xyz.z);
	a = tmp;
}




/*-------------------------------------


			OPENGL UTILITY


---------------------------------------*/

const cmat4  IDENTITY_MAT = {{1,0,0,0,
0,1,0,0,     0,0,1,0,     0,0,0,1}};

float  cotangent(float ang){
	return 1.0/tan(ang);
}

float  deg2rad(float deg){
	return deg*cPI/180.0;
}

float  rad2deg(float rad){
	return rad*cPI/180.0;
}

void  RotateAboutX(cmat4& m,float ang){
	cmat4 rotation = IDENTITY_MAT;
	float sine = sin(ang);
	float cosine = cos(ang);

	rotation[1][1] = cosine;
	rotation[1][2] = -sine;
	rotation[2][1] = sine;
	rotation[2][2] = cosine;

	m = m * rotation;
}

void  RotateAboutY(cmat4& m, float ang) {
	cmat4 rotation = IDENTITY_MAT;
	float sine = sin(ang);
	float cosine = cos(ang);

	rotation[0][0] = cosine;
	rotation[0][2] = -sine;
	rotation[2][0] = sine;
	rotation[2][2] = cosine;

	m = m * rotation;
}

void  RotateAboutZ(cmat4& m, float ang) {
	cmat4 rotation = IDENTITY_MAT;
	float sine = sin(ang);
	float cosine = cos(ang);

	rotation[0][0] = cosine;
	rotation[0][1] = -sine;
	rotation[1][0] = sine;
	rotation[1][1] = cosine;

	m = m * rotation;
}

void  ScaleMatrix(cmat4& m, cfloat3 x){
	cmat4 scale = IDENTITY_MAT;
	scale[0][0] = x.x;
	scale[1][1] = x.y;
	scale[2][2] = x.z;
	m = m*scale;
}

void  TranslateMatrix(cmat4& m, cfloat3 x){
	cmat4 trans = IDENTITY_MAT;
	trans[3][0] = x.x;
	trans[3][1] = x.y;
	trans[3][2] = x.z;
	m = m*trans;
}

cmat4 CreateProjectionMatrix(float fovy, float aspect_ratio, float near_plane, float far_plane){
	cmat4 out = {{0}};
	const float 
		y_scale = cotangent(deg2rad(fovy/2)),
		x_scale = y_scale / aspect_ratio,
		frustum_length = far_plane - near_plane;
	out[0][0] = x_scale;
	out[1][1] = y_scale;
	out[2][2] = -((far_plane + near_plane)/frustum_length);
	out[2][3] = -1;
	out[3][2] = -((2*near_plane * far_plane)/frustum_length);
	
	return out;
}

inline cfloat3  cross(cfloat3& a, cfloat3& b) {
	cfloat3 tmp;
	tmp.x = a.y*b.z - a.z*b.y;
	tmp.y = a.z*b.x - a.x*b.z;
	tmp.z = a.x*b.y - a.y*b.x;
	return tmp;
}




float  angle(cfloat3& a,cfloat3& b){
	return acos( dot(a,b) / sqrt(dot(a,a)) /sqrt(dot(b,b)) );
}

cmat4  cCamera::CameraRotateMat(){

	cfloat3 right = cross(forceup,dir);
	right = right * (1/sqrt(dot(right,right)));
	this->right = right;
	up = cross(dir,right); //update up

	cmat4 tmp = IDENTITY_MAT;
	tmp[0][0] = right.x;
	tmp[1][0] = right.y;
	tmp[2][0] = right.z;
	tmp[3][0] = 0; //new x

	tmp[0][1] = up.x;
	tmp[1][1] = up.y;
	tmp[2][1] = up.z;
	tmp[3][1] = 0; //new y

	tmp[0][2] = dir.x;
	tmp[1][2] = dir.y;
	tmp[2][2] = dir.z;
	tmp[3][2] = 0; //new z

	tmp[0][3] = 0;
	tmp[1][3] = 0;
	tmp[2][3] = 0;
	tmp[3][3] = 1; //new w
	return tmp;
}

void  cCamera::CameraTranslateMat(cmat4& rotation,cfloat3& pos){
	cmat4 tmp = IDENTITY_MAT;
	tmp[3][0]=-pos.x;
	tmp[3][1]=-pos.y;
	tmp[3][2]=-pos.z;

	rotation =  tmp * rotation;
}

void cCamera::AdvanceCamera(float dt){
	if (wasd.lenth==0) //no input
	{
		if (velf > 0.001)
			velf -= velf*10*dt;
	}
	else{
		velf += (velMax - velf) * 5 * dt;
		getVelDir();
	}

	if (velf < 0.001 && dot(rotatexy, rotatexy)<0.001 && dot(transxy,transxy)<0.001) {
		return;
	}
	cfloat3 dx = vel * velf * dt;

	//translate
	pos = pos+dx;
	target = target+dx;

	//rotate
	
	rotatexy  = rotatexy * 0.01;
	target = target - up*rotatexy.y;
	target = target + right*rotatexy.x;

	transxy = transxy * 0.01;
	dx = up * (transxy.y) + right*(-transxy.x);
	dx = dx*100;
	pos = pos+dx;
	target = target+dx;

	lookat(pos, target);
}