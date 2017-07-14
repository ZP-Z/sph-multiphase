
#include <assert.h>
#include <stdio.h>
#include <conio.h>

#include "fluid_system.h"
#include "fluid_system_host.cuh"

double scaleP,scaleP3,scaledis;

cfloat3 volumes[10];

int example = 2; 

float unitMatrix[9] = {1,0,0,    0,1,0,     0,0,1};

FluidSystem::FluidSystem ()
{
	mPos = 0x0;
	mClr = 0x0;
	mIsBound = 0x0;
	mVel = 0x0;
	mVelEval = 0x0;
	mPressure = 0x0;
	mDensity = NULL;
	mForce = NULL;

	m_alpha = NULL;
	m_alpha_pre = NULL;
	m_pressure_modify = NULL;
	m_vel_phrel = NULL;

	m_restMass = NULL;
	m_restDensity = NULL;
	m_visc = NULL;
	m_velxcor = NULL;
	m_alphagrad = NULL;

}


void FluidSystem::Setup ( bool bStart )
{
	
	ResetParameters();

	ParseXML();

    //SetupMpmCase();
	SetupSimpleSphCase();

	SetupGridAllocate ();	// Setup grid

	SetupDevice();
	
	//BeforeFirstRun();
}

void FluidSystem::Exit ()
{
	free ( mPos );
	free ( mClr );
	free (mIsBound);
	free ( mVel );
	free ( mVelEval );
	free ( mPressure );
	free ( mDensity );
	free ( mForce );
	
	//multi fluid
	free (m_alpha);
	free (m_alpha_pre);
	free (m_pressure_modify);
	free (m_vel_phrel);
	free (m_restMass);
	free (m_restDensity);
	free (m_visc);
	free (m_velxcor);
	free (m_alphagrad);
	free (MF_type);
	free (MF_tensor);

	FluidClearCUDA();

	cudaExit (0,0);
}




// Setting Parameters

void FluidSystem::ResetParameters()
{
	frameNo = 0;
	outputNo = 0;
	pointNum = 0;
	lastTime = 0.0f;
}

void FluidSystem::AllocateParticles(int cnt)
{
	int nump = 0;		// number to copy from previous data

	displayBuffer = new displayPack[EMIT_BUF_RATIO * cnt];
	calculationBuffer = new calculationPack[EMIT_BUF_RATIO*cnt];

	mPos = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*sizeof(cfloat3));
	mClr = (DWORD*)malloc(EMIT_BUF_RATIO*cnt*sizeof(DWORD));
	mColor = new cfloat4[EMIT_BUF_RATIO*cnt];

	mIsBound = (int*)malloc(cnt*sizeof(int));
	mVel = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*sizeof(cfloat3));
	mVelEval = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*sizeof(cfloat3));
	mPressure = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));
	mDensity = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));
	mForce = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*sizeof(cfloat3));

	//multi fluid

	m_alpha = (float*)malloc(EMIT_BUF_RATIO*cnt*MAX_FLUIDNUM*sizeof(float));
	m_alpha_pre = (float*)malloc(EMIT_BUF_RATIO*cnt*MAX_FLUIDNUM*sizeof(float));
	m_pressure_modify = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));
	m_vel_phrel = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*MAX_FLUIDNUM*sizeof(cfloat3));

	m_restMass = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));
	m_restDensity = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));
	m_visc = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));

	m_velxcor = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*sizeof(cfloat3));
	m_alphagrad = (cfloat3*)malloc(EMIT_BUF_RATIO*cnt*MAX_FLUIDNUM*sizeof(cfloat3));

	//Project U

	MF_type = (int*)malloc(EMIT_BUF_RATIO*cnt*sizeof(int));
	MF_tensor = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(int)*9);

	MF_id = (int*)malloc(EMIT_BUF_RATIO*cnt*sizeof(int));
	MF_idTable = (int*)malloc(EMIT_BUF_RATIO*cnt*sizeof(int));
	MF_pepsilon = (float*)malloc(EMIT_BUF_RATIO*cnt*sizeof(float));

	//End Project U
	printf("particle buffer %d\n", cnt);
}

void FluidSystem::SetupSpacing ()
{
	// Determine spacing from density
	pRealDist = pow( hostCarrier.mass/ hostCarrier.restdensity, 1/3.0);
	pSpacing = pRealDist /  hostCarrier.simscale;
	
	printf ( "Add Particles. Mass: %f, Density: %f, Spacing: %f PDist: %f\n", 
		hostCarrier.mass,
		hostCarrier.restdensity,
		pSpacing,
		pRealDist);

}

void FluidSystem::SetupGridAllocate ()
{
	float worldcellsize = hostCarrier.cellsize / hostCarrier.simscale;
	
	hostCarrier.gridmin = hostCarrier.volmin;
	hostCarrier.gridmax = hostCarrier.volmax;
	hostCarrier.gridsize = hostCarrier.gridmax - hostCarrier.gridmin;
	hostCarrier.gridres.x = ceil(hostCarrier.gridsize.x / worldcellsize);
	hostCarrier.gridres.y = ceil(hostCarrier.gridsize.y / worldcellsize);
	hostCarrier.gridres.z = ceil(hostCarrier.gridsize.z / worldcellsize);
	hostCarrier.gridsize.x = hostCarrier.gridres.x * worldcellsize;
	hostCarrier.gridsize.y = hostCarrier.gridres.y * worldcellsize;
	hostCarrier.gridsize.z = hostCarrier.gridres.z * worldcellsize;
	hostCarrier.gridIdfac.x =  1.0 / hostCarrier.gridsize.x * hostCarrier.gridres.x; 
	hostCarrier.gridIdfac.y =  1.0 / hostCarrier.gridsize.y * hostCarrier.gridres.y;
	hostCarrier.gridIdfac.z =  1.0 / hostCarrier.gridsize.z * hostCarrier.gridres.z;

	hostCarrier.gridtotal = hostCarrier.gridres.x * hostCarrier.gridres.y * hostCarrier.gridres.z;
	hostCarrier.searchnum = 3;
	hostCarrier.neighbornum = 27;
}

void FluidSystem::SetupDevice() {

	CarryParam(hostCarrier);

	FluidSetupCUDA(hostCarrier); //allocate buffers
	FluidParamCUDA(hostCarrier); //some parameters

	TransferToCUDA();		// Initial transfer

}

void FluidSystem::BeforeFirstRun() {
	//Sorting
	InitialSortCUDA(0x0, 0x0, 0x0);
	SortGridCUDA(0x0);
	CountingSortFullCUDA_(0x0);
	initSPH(m_restDensity, MF_type);


	InitSolid(); //clear tensor buffer
}



//Adding Particles


//void FluidSystem::setupSPHexample(){
//	ParseXML_Bound("BoundInfo",1);
//	example = 20;
//	double particleVisc =  m_Param[PVISC];
//	
//	//BI2Reader* bi2readers[10]; //ten pointers to build bi2reader dynamically, in use now
//	//char biname[200];
//	//sprintf(biname,".\\BoundaryOut\\BoundaryOut\\Boundary.bi2");
//	//bi2readers[0] = new BI2Reader(biname);
//	//bi2readers[0]->GetInfo(false);
//	//bi2readers[0]->PrintInfo();
//	
//	//parse the xml and adjust some parameters according to scaleP
//	ParseMFXML ( "MultiScene", example, true );
//	coK = coE/3./(1-2*coV);
//	coG = coE / 2. / (1 + coV);
//	solid_coK = solid_coG / 3. / (1 - 2 * solid_coV);
//	solid_coG = solid_coG / 2. / (1 + solid_coV);
//
//	//adjust the parametres according to the scale parameter
//	scaleP3 = pow(scaleP,1.0/3.0);
//	m_Param[PMASS]/=scaleP;
//	m_Param[PBMASS]/=scaleP;
//	m_Param[PSMOOTHRADIUS]/=scaleP3;
//	m_Param[PRADIUS]/=scaleP3;
//	m_Param [PNUM]*=scaleP;
//	//m_Param[PNUM] += bi2readers[0]->info.nbound;
//
//	m_Param [PGRIDSIZE] = 2*m_Param[PSMOOTHRADIUS] / m_Param[PGRID_DENSITY];
//	
//	for(int k=0; k<3; k++){
//		m_fluidPMass[k]   = m_Param[PMASS] * _massRatio[k];
//		m_fluidDensity[k] = 600.0 * _densityRatio[k];
//		m_fluidVisc[k]    = particleVisc * _viscRatio[k];
//	}
//
//	AllocateParticles ( m_Param[PNUM] );
//	SetupKernels ();
//	SetupSpacing ();
//
//	m_Vec [ PINITMIN ].Set (volumes[0].x,volumes[0].y,volumes[0].z);
//	m_Vec [ PINITMAX ].Set (volumes[1].x,volumes[1].y,volumes[1].z);
//	//SetupMfAddVolume ( m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING], cfloat3(0,0.1,0),1);
//	SetupMfAddDeformVolume ( m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING]*spacing_multipler, cfloat3(0,0.1,0),0);
//	
//	int open = 0;
//	if(open){
//		m_Vec [ PINITMIN ].Set (volumes[2].x,volumes[2].y,volumes[2].z);
//		m_Vec[PINITMAX].Set(volumes[3].x, volumes[3].y, volumes[3].z);
//		//SetupMfAddDeformVolume(m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING] * solid_spacing_multipler, cfloat3(0, 0.1, 0), 1);
//		SetupMfAddDeformVolume ( m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING] * spacing_multipler, cfloat3(0, 0.1, 0), 1);
//
//		m_Vec [ PINITMIN ].Set (volumes[4].x,volumes[4].y,volumes[4].z);
//		m_Vec [ PINITMAX ].Set (volumes[5].x,volumes[5].y,volumes[5].z);
//		SetupMfAddVolume ( m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING]*1.1, cfloat3(0,0.1,0),2);
//
//		m_Vec [ PINITMIN ].Set (volumes[6].x,volumes[6].y,volumes[6].z);
//		m_Vec [ PINITMAX ].Set (volumes[7].x,volumes[7].y,volumes[7].z);
//		SetupMfAddVolume ( m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING]*1.1, cfloat3(0,0.1,0),2);
//
//		m_Vec [ PINITMIN ].Set (volumes[8].x,volumes[8].y,volumes[8].z);
//		m_Vec [ PINITMAX ].Set (volumes[9].x,volumes[9].y,volumes[9].z);
//		SetupMfAddVolume ( m_Vec[PINITMIN], m_Vec[PINITMAX], m_Param[PSPACING]*1.1, cfloat3(0,0.1,0),2);
//	}
//	//SetupAddBound(*bi2readers[0],1);
//	m_maxAllowedPoints = mMaxPoints;
//
//	//Emit Parameters
//	m_Vec[PEMIT_RATE] = cfloat3(75,10000*scaleP,36);
//	m_Vec[PEMIT_SPREAD] = cfloat3(0,0,0);
//	m_Vec[PEMIT_ANG] = cfloat3(0,90,0.6);
//	m_Vec[PEMIT_POS] = cfloat3(0,60,-17);
//}
//
//void FluidSystem::SetupMpmCase(){
//	ParseXML_Bound("BoundInfo",1);
//	example = 20;
//	double particleVisc =  m_Param[PVISC];
//	
//	//BI2Reader* bi2readers[10]; //ten pointers to build bi2reader dynamically, in use now
//	//char biname[200];
//	//sprintf(biname,".\\BoundaryOut\\BoundaryOut\\Boundary.bi2");
//	//bi2readers[0] = new BI2Reader(biname);
//	//bi2readers[0]->GetInfo(false);
//	//bi2readers[0]->PrintInfo();
//	
//	//parse the xml and adjust some parameters according to scaleP
//	ParseMFXML ( "MultiScene", example, true );
//	coK = coE/3./(1-2*coV);
//	coG = coE / 2. / (1 + coV);
//	solid_coK = solid_coG / 3. / (1 - 2 * solid_coV); //coG=coE
//	solid_coG = solid_coG / 2. / (1 + solid_coV);     //coG=coG
//
//	//adjust the parametres according to the scale parameter
//	scaleP3 = pow(scaleP,1.0/3.0);
//	m_Param[PMASS]/=scaleP;
//	m_Param[PBMASS]/=scaleP;
//	m_Param[PSMOOTHRADIUS]/=scaleP3;
//	m_Param[PRADIUS]/=scaleP3;
//	m_Param [PNUM]*=scaleP;
//	//m_Param[PNUM] += bi2readers[0]->info.nbound;
//
//	m_Param [PGRIDSIZE] = 2*m_Param[PSMOOTHRADIUS] / m_Param[PGRID_DENSITY];
//	
//	for(int k=0; k<3; k++){
//		m_fluidPMass[k]   = m_Param[PMASS] * _massRatio[k];
//		m_fluidDensity[k] = 600.0 * _densityRatio[k];
//		m_fluidVisc[k]    = particleVisc * _viscRatio[k];
//	}
//	
//	AllocateParticles ( m_Param[PNUM] );
//    mMaxPoints = m_Param[PNUM];
//    m_maxAllowedPoints = mMaxPoints;
//
//	SetupKernels ();
//	SetupSpacing ();
//
//	//SetupMfAddVolume(volumes[0], volumes[1], m_Param[PSPACING], cfloat3(0, 0.1, 0), 0);
//    SetupMfAddVolume(volumes[2], volumes[3], m_Param[PSPACING]*1.1, cfloat3(0, 0.1, 0), 1);
//
//	//SetupMfAddSolid(volumes[0], volumes[1], m_Param[PSPACING], cfloat3(0, 0.1, 0), 0);
//	//SetupMfAddSolid(volumes[2], volumes[3], m_Param[PSPACING], cfloat3(0, 0.1, 0), 0);
//
//	//SetupAddBound(*bi2readers[0],1);
//	
//}

void FluidSystem::SetupSimpleSphCase(){
	
	//adjust the parametres according to the scale parameter
	scaleP3 = pow(scaleP,1.0/3.0);
	hostCarrier.mass /= scaleP;
	hostCarrier.smoothradius /= scaleP3;
	hostCarrier.radius /= scaleP3;
	hostCarrier.num *= scaleP;

	hostCarrier.cellsize = hostCarrier.smoothradius;
	
	for(int k=0; k<fluidnum; k++){
		hostCarrier.massArr[k] = hostCarrier.mass * massratio[k];
		hostCarrier.densArr[k] = hostCarrier.restdensity * densratio[k];
		hostCarrier.viscArr[k] = hostCarrier.viscosity * viscratio[k];
	}

	AllocateParticles ( maxPointNum );
	SetupSpacing ();
	
	SetupMfAddVolume(fvmin[0], fvmax[0], pSpacing, cfloat3(0, 0, 0), 2);
	//SetupAddSolid(volumes[0], volumes[1], m_Param[PSPACING], 0);
    //SetupAddBubble(volumes[2], volumes[3], m_Param[PSPACING]*1.2, 1);

    //LoadBoundary("cfg\\boundary.cfg");
	//LoadBoundary("cfg\\Cup.cfg");

	hostCarrier.num = pointNum;
}

//Solver
//
//void FluidSystem::RunSimulateMultiCUDAFull ()
//{
//	cTime start;
//	
//	//start.SetSystemTime ( ACC_NSEC );
//
//	//------------------Sorting-------------------
//	InitialSortCUDA( 0x0, 0x0, 0x0 );			
//	SortGridCUDA( 0x0 );
//	CountingSortFullCUDA_( 0x0 );
//	record ( PTIME_SORT, "Full Sort CUDA", start );
//
//	start.update();
//
//	//-------------Compute Pressure--------------
//	MfComputePressureCUDA();                                          //case 3,5
//	record ( PTIME_PRESS, "Compute Pressure CUDA", start);
//	start.update();
//
//	//----------Compute Drift Velocity-----------
//	MfComputeDriftVelCUDA();                                          //case 1-diff
//	record ( PTIMEDRIFTVEL, "Drift Velocity CUDA", start);
//	start.update();
//
//	//------------Compute Alpha Change--------------
//	MfComputeAlphaAdvanceCUDA();									//case 1
//	record ( PTIMEALPHA, "Alpha Advance CUDA", start );
//	start.update();
//
//	//---------Compute Pressure Correction-------------
//	MfComputeCorrectionCUDA();                                        //case5
//	record ( PTIMECORR, "Alpha Correction and Pressure CUDA", start );		
//	start.update();
//
//	//-----------------Compute Force-------------------
//	//MfComputeForceCUDA ();                                            //case 5
//	ComputeForceCUDA_ProjectU();
//	record ( PTIME_FORCE, "Force CUDA", start );
//	start.update();
//
//	//---------------Particle Advance---------------------
//	MfAdvanceCUDA ( m_Time, m_DT, m_Param[PSIMSCALE] );	             //case 3 5
//	record ( PTIME_ADVANCE, "Advance CUDA", start );
//	
//	TransferFromCUDA ();	// return for rendering
//}
//
//void FluidSystem::RunMpmSolution(){
//    
//    //mint::Time start;
//	//start.SetSystemTime ( ACC_NSEC );
//
//	InitialSortCUDA( 0x0, 0x0, 0x0 );			
//	SortGridCUDA( 0x0 );
//	CountingSortFullCUDA_( 0x0 );
//	//record ( PTIME_SORT, "Full Sort CUDA", start );
//
//	MfComputePressureCUDA();
//	//record ( PTIME_PRESS, "Compute Pressure CUDA", start );
//	//start.SetSystemTime ( ACC_NSEC );
//
//	MfComputeDriftVelCUDA();
//	//record ( PTIMEDRIFTVEL, "Drift Velocity CUDA", start );
//	//start.SetSystemTime ( ACC_NSEC );
//
//	MfComputeAlphaAdvanceCUDA();
//	//record ( PTIMEALPHA, "Alpha Advance CUDA", start );
//	//start.SetSystemTime ( ACC_NSEC );
//
//	MfComputeCorrectionCUDA();
//	//record ( PTIMECORR, "Alpha Correction and Pressure CUDA", start );		
//	//start.SetSystemTime ( ACC_NSEC );
//
//	//ComputeForceCUDA_ProjectU();
//    ComputeMpmForce();                                                       // <<<<<<<<<<<<<<  MPM Function
//
//    //record ( PTIME_FORCE, "Force CUDA", start );
//	//start.SetSystemTime ( ACC_NSEC );
//
// //   MpmParticleCollision();
//	//MfAdvanceCUDA ( m_Time, m_DT, m_Param[PSIMSCALE] );
//	//record ( PTIME_ADVANCE, "Advance CUDA", start );
//	//
//	TransferFromCUDA ();
//}




//
//void FluidSystem::RunSolid(){
//	cTime start;
//	
//	//------------------Sorting-------------------
//	InitialSortCUDA(0x0, 0x0, 0x0);
//	SortGridCUDA(0x0);
//	CountingSortFullCUDA_(0x0);
//	record(PTIME_SORT, "Full Sort CUDA", start);
//	start.update();
//	//Velocity Gradient, Strain, Stress
//	ComputeSolidTensor();
//	record(PTIME_PRESS, "Stress", start);
//	start.update();
//	//Force
//	ComputeSolidForce();
//	record(PTIME_FORCE, "Force", start);
//	start.update();
//	//Advance
//	MfAdvanceCUDA(m_Time, m_DT, m_Param[PSIMSCALE]);
//	record(PTIME_ADVANCE, "Advance CUDA", start);
//
//	TransferFromCUDA();
//}


cTime start;

void FluidSystem::ClearTimer(){
	timerStrs.clear();
	timerVals.clear();
	start.update();
}

void FluidSystem::CheckTimer(const char* msg){
	timerVals.push_back(cTime() - start);
	timerStrs.push_back(msg);
	start.update();
}

void FluidSystem::Run (int width, int height)
{
	
	runMode = 1;

	switch( runMode){
	case 0:	
		//RunMpmSolution();
		break;
	case 1:
		RunSimpleSPH();
		break;
	case 2:
		//RunSolid();
		break;
	}

	if ( bSnapshot && frameNo % 20 ==0){
		CaptureVideo ( width, height );
	}

	//if ( frameNo % framespacing ==0 ){
	//	outputFile();
	//}

	frameNo++;
}

extern bufList	fbuf;
void FluidSystem::RunSimpleSPH() {

	ClearTimer();

	//------------------Sorting-------------------
	InitialSortCUDA(0x0, 0x0, 0x0);
	SortGridCUDA(0x0);
	CountingSortFullCUDA_(0x0);
	CheckTimer("sorting");

	//-------------Compute Pressure--------------
	MfComputePressureCUDA();
	CheckTimer("pressure");

	//-----------------Compute Force-------------------
	MfComputeForceCUDA();
	CheckTimer("force");
	//---------------Particle Advance---------------------
	MfAdvanceCUDA(0, hostCarrier.dt, hostCarrier.simscale);
	CheckTimer("advance");

	TransferFromCUDA();	// return for rendering
	TransferFromCUDA(fbuf);

	for (int i=0; i<pointNum; i++) {
		if (MF_id[i]==1000)
			mColor[i].Set(1, 0, 0, 0.9);
		else
			mColor[i].w = 0.5;
	}
}




void FluidSystem::EmitParticles (int cat)
{
	int emitInterval = 20;
	int oldPointNum = pointNum;
    if (pointNum == maxPointNum)
        return;
	std::vector<cfloat3> plist;

	if (  frameNo % emitInterval == 0 ) {
		ParticleEmitter emitter;
		emitter.dir.Set(0,-1,0);
		emitter.position.Set(0,100,0);
		emitter.vel = 5;
		emitter.radius = 10;
		emitter.spacing = pSpacing;
		emitter.EmitImmediate(plist);

		InsertParticles(plist, cat);
	}

	EmitUpdateCUDA(oldPointNum,pointNum);
}


void ParticleEmitter::EmitImmediate(std::vector<cfloat3>& plist){
	//plist.clear();
	cfloat3 start;
	
	cfloat3 inletCenter;
	cmat3 rotationMat;
	rotationMat.Set(unitMatrix);
	cfloat3 newx = cross(dir, cfloat3(0,1,0));
	if( dot(newx,newx)> 0.001) //not aligned
	{
		newx = newx / sqrt(dot(newx,newx));
		cfloat3 newz = cross(newx,dir);
		printf("%f\n", dot(newz,newz));//should be 1
		rotationMat[0][0] = newx.x; rotationMat[0][1] = dir.x; rotationMat[0][2] = newz.x;
		rotationMat[1][0] = newx.y; rotationMat[1][1] = dir.y; rotationMat[1][2] = newz.y;
		rotationMat[2][0] = newx.z; rotationMat[2][1] = dir.z; rotationMat[2][2] = newz.z;
	}

	int res = radius / spacing;
	for(int i=-res; i<=res; i++){
		for(int j=-res;j<=res;j++){
			cfloat3 tmp( i*spacing, 0, j*spacing);
			if( dot(tmp,tmp)>radius*radius)
				continue;
			mvprod(rotationMat, tmp, tmp);
			tmp += position;
			plist.push_back(tmp);
		}
	}
}
void FluidSystem::InsertParticles(std::vector<cfloat3> plist, int cat){
	for(int i=0; i<plist.size();i++){
		int p = AddParticle();

	}
}

XMLElement* pele=NULL;
inline float XMLGetFloat(const char* name){
	float tmp=0;
	if(pele!=NULL){
		XMLElement* attr = pele->FirstChildElement(name);
		if(attr)
			attr->QueryFloatText(&tmp);
	}
	return tmp;
}

inline int XMLGetInt(const char* name) {
	int tmp=0;
	if (pele!=NULL) {
		pele->FirstChildElement(name)->QueryIntText(&tmp);
	}
	return tmp;
}

inline cfloat3& XMLGetFloat3(const char* name) {
	if (pele!=NULL) {
		return QueryFloat3(pele->FirstChildElement(name));
	}else
		return cfloat3(0,0,0);
}


#include <string>
void XMLGetFloatN(float* buffer, int size, const char* name) {
	if (pele==NULL)
		return;
	const char* str = pele->FirstChildElement(name)->GetText();
	std::string fmstr = "";
	for(int i=0;i<size-1;i++)
		fmstr = fmstr+"%f,";
	fmstr = fmstr+"%f";
	sscanf(str,fmstr.c_str(),&buffer[0], &buffer[1], &buffer[2]);
}

void FluidSystem::ParseXML(){
	tinyxml2::XMLDocument doc;
	int result = doc.LoadFile("scene.xml");
	printf("return by open xml %d\n",result);

	XMLElement* fluidElement = doc.FirstChildElement("Fluid");
	XMLElement* sceneElement = doc.FirstChildElement("MultiScene");
	
	if(!fluidElement || !sceneElement)
	{
		printf("missing fluid/scene xml node");
		return;
	}

	int tmp;
	int caseid=2;
	while(true){
		sceneElement -> QueryIntAttribute("id",&tmp);
		if(tmp == caseid)
			break;
		else if(sceneElement->NextSiblingElement() != NULL){
			sceneElement = sceneElement->NextSiblingElement();
		}
		else{
			printf("scene not found.\n");
			return;
		}
	}

	//some general parameter
	pele = fluidElement;
	hostCarrier.intstiff = XMLGetFloat("IntStiff");
	hostCarrier.extstiff = XMLGetFloat("ExtStiff");
	hostCarrier.extdamp = XMLGetFloat("ExtDamp");
	hostCarrier.acclimit = XMLGetFloat("AccelLimit");
	hostCarrier.vlimit = XMLGetFloat("VelLimit");
	hostCarrier.gravity = XMLGetFloat3("Gravity");
	hostCarrier.volmin = XMLGetFloat3("VolMin");
	hostCarrier.volmax = XMLGetFloat3("VolMax");

	pele = sceneElement;
	//scene specified
	maxPointNum = XMLGetInt("Num");
	scaleP = XMLGetFloat("ScaleP");
	fluidnum = XMLGetInt("FluidCount");
	fvmin[0] = XMLGetFloat3("VolMin0");
	fvmax[0] = XMLGetFloat3("VolMax0");
	XMLGetFloatN(massratio, 3, "MassRatio");
	XMLGetFloatN(densratio,3,"DensityRatio");
	XMLGetFloatN(viscratio, 3, "ViscRatio");

	hostCarrier.simscale = XMLGetFloat("SimScale");
	hostCarrier.dt = XMLGetFloat("DT");
	hostCarrier.smoothradius = XMLGetFloat("SmoothRadius");
	hostCarrier.viscosity = XMLGetFloat("Viscosity");
	hostCarrier.mass = XMLGetFloat("Mass");
	hostCarrier.restdensity = XMLGetFloat("RestDensity");
	hostCarrier.radius = XMLGetFloat("Radius");
	hostCarrier.extstiff = XMLGetFloat("ExtStiff");
	hostCarrier.intstiff = XMLGetFloat("IntStiff");
	hostCarrier.extdamp = XMLGetFloat("ExtDamp");
	hostCarrier.softminx = XMLGetFloat3("SoftBoundMin");
	hostCarrier.softmaxx = XMLGetFloat3("SoftBoundMax");
}

void FluidSystem::TransferToCUDA(bufList& fbuf) {
	cudaMemcpy(fbuf.displayBuffer,displayBuffer,sizeof(displayPack)*pointNum, cudaMemcpyHostToDevice);
	//cudaMemcpy(fbuf.calcBuffer, calculationBuffer, sizeof(calculationPack)*pointNum, cudaMemcpyHostToDevice);
	//cudaMemcpy(fbuf.intmBuffer, intmBuffer, sizeof(IntermediatePack)*pointNum, cudaMemcpyHostToDevice);
}

void FluidSystem::TransferFromCUDA(bufList& fbuf) {
	cudaMemcpy( displayBuffer, fbuf.displayBuffer, sizeof(displayPack)*pointNum, cudaMemcpyDeviceToHost);
	//cudaMemcpy( calculationBuffer, fbuf.calcBuffer, sizeof(calculationPack)*pointNum, cudaMemcpyDeviceToHost);
	//cudaMemcpy(intmBuffer, fbuf.intmBuffer, sizeof(IntermediatePack)*pointNum, cudaMemcpyDeviceToHost);
}