
#include <assert.h>
#include <stdio.h>
#include <conio.h>

#include "fluid_system.h"
#include "fluid_system_host.cuh"
#include "MpmSolver.cuh"

double scaleP,scaleP3,scaledis;

cfloat3 volumes[10];

int example = 2; 

float unitMatrix[9] = {1,0,0,    0,1,0,     0,0,1};

extern bufList	fbuf;
extern FluidParams		fcuda;
ParamCarrier hostCarrier;



FluidSystem::FluidSystem ()
{
}

#define RUN_SPH 1
#define RUN_SOLID 2
#define RUN_IISPH 3
#define RUN_MPM 4

#define SPH_CASE 1
#define SOLID_CASE 2

void FluidSystem::Setup ()
{
	
	ResetParameters();
	
	runMode = RUN_MPM;
	glutSetWindowTitle("Fluid SPH");

	int setupCase = SPH_CASE;
	switch (setupCase) {
	case SPH_CASE:
		SetupSimpleSphCase();
		break;
	case SOLID_CASE:
		SetupSolidCase();
		break;
	}

	SetupGrid();
	SetupDevice();
	
	BeforeFirstStep();
	SetupMPMGrid();
	//InitializeSolid_CUDA();
}

void FluidSystem::BeforeFirstStep() {
	
	GetParticleIndexCUDA();
	GetGridListCUDA();
	RearrageDataCUDA();

	printf("Initialize Boundary Density Number.\n");
	ComputeBoundaryDensity();
}

void FluidSystem::Exit ()
{
	delete displayBuffer;
	delete calculationBuffer;

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
	bSnapshot = false;
}

void FluidSystem::AllocateParticles(int cnt)
{
	int nump = 0;		// number to copy from previous data

	displayBuffer = new displayPack[EMIT_BUF_RATIO * cnt];
	calculationBuffer = new calculationPack[EMIT_BUF_RATIO*cnt];

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

void FluidSystem::SetupGrid ()
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

	for(int i=0; i<3; i++)
		for(int j=0; j<3; j++)
			for(int k=0; k<3; k++)
				hostCarrier.neighborid[i*9+j*3+k]=(i-1)*hostCarrier.gridres.x*hostCarrier.gridres.z+(j-1)*hostCarrier.gridres.x+(k-1);

}

void FluidSystem::SetupDevice() {

	FluidSetupCUDA(hostCarrier); //allocate buffers
	//FluidParamCUDA(hostCarrier); //some parameters
	
	CarryParam(hostCarrier);
	TransferToCUDA(fbuf);

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

void FluidSystem::SetSolidModulus(float E, float v) {
	solidE = E;
	solidv = v;
	solidG = solidE / 2/ (1+solidv);
	solidK = solidE / 3/ (1-2*solidv);
	printf("Setting G:%f K:%f\n",solidG,solidK);
}



void FluidSystem::SetupSimpleSphCase(){
	
	ParseXML(2);

	//adjust the parametres according to the scale parameter
	/*scaleP3 = pow(scaleP,1.0/3.0);
	hostCarrier.mass /= scaleP;
	hostCarrier.smoothradius /= scaleP3;
	hostCarrier.radius /= scaleP3;
	hostCarrier.num *= scaleP;*/

	hostCarrier.cellsize = hostCarrier.smoothradius*2;
	hostCarrier.kpoly6 = 315.0f / (64.0f * 3.141592 * pow(hostCarrier.smoothradius, 9.0f));
	hostCarrier.kspiky =  15 / (3.141592 * pow(hostCarrier.smoothradius, 6.0f));
	hostCarrier.kspikydiff = -45.0f / (3.141592 * pow(hostCarrier.smoothradius, 6.0f));
	hostCarrier.klaplacian = 45.0f / (3.141592 * pow(hostCarrier.smoothradius, 6.0f));
	hostCarrier.kspline = 1.0f/3.141593f/pow(hostCarrier.smoothradius, 3.0f);

	for(int k=0; k<fluidnum; k++){
		hostCarrier.massArr[k] = hostCarrier.mass * massratio[k];
		hostCarrier.densArr[k] = hostCarrier.restdensity * densratio[k];
		hostCarrier.viscArr[k] = hostCarrier.viscosity * viscratio[k];
	}

	AllocateParticles ( maxPointNum );
	SetupSpacing ();
	
	
	SetSolidModulus(10000000000, 0.3);
	AddGranularVolume(fvmin[0], fvmax[0], pSpacing, cfloat3(0, 0, 0));
	//AddFluidVolume(fvmin[0], fvmax[0], pSpacing, cfloat3(0, 0, 0),2);
	SetSolidModulus(10000000000, 0.3);
	AddFluidVolume(fvmin[1], fvmax[1], pSpacing, cfloat3(5, 15, 0), 2);
	
	if (false) {
		SetSolidModulus(500000000, 0.3);
		BI2Reader* bi2readers;
		bi2readers = new BI2Reader("boundaryout\\bunny\\bunny_def.bi2");
		bi2readers->GetInfo(false);
		bi2readers->PrintInfo();
		AddShape(*bi2readers, TYPE_ELASTIC);
		delete bi2readers;
	}
	

    LoadBoundary("cfg\\boundary-wide.cfg");
	
	hostCarrier.num = pointNum;
}




void FluidSystem::SetupSolidCase() {

	ParseXML(2);

	hostCarrier.cellsize = hostCarrier.smoothradius;

	for (int k=0; k<fluidnum; k++) {
		hostCarrier.massArr[k] = hostCarrier.mass * massratio[k];
		hostCarrier.densArr[k] = hostCarrier.restdensity * densratio[k];
		hostCarrier.viscArr[k] = hostCarrier.viscosity * viscratio[k];
	}

	AllocateParticles(maxPointNum);
	SetupSpacing();

	AddDeformableVolume(fvmin[0], fvmax[0], pSpacing, 2);
	
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




//timer for simulation stages
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





void FluidSystem::Run()
{
	
	

	switch( runMode){
	case RUN_SPH:
		RunSimpleSPH();
		break;
	case RUN_SOLID:
		RunSolid();
		break;
	case RUN_IISPH:
		RunIISPH();
		break;
	case RUN_MPM:
		RunMPM();
		break;
	}

	//if (frameNo==20) {
	//EmitParticles(2);
	//}

	//if ( frameNo % framespacing ==0 ){
	//	outputFile();
	//}

	frameNo++;
}


void FluidSystem::RunSimpleSPH() {

	ClearTimer();

	//------------------Sorting-------------------
	GetParticleIndexCUDA();
	GetGridListCUDA();
	RearrageDataCUDA();
	CheckTimer("sorting");

	//-------------Compute Pressure--------------
	MfComputePressureCUDA();
	CheckTimer("pressure");

	//-----------------Compute Force-------------------
	MfComputeForceCUDA();
	CheckTimer("force");
	//---------------Particle Advance---------------------
	MfAdvanceCUDA();
	CheckTimer("advance");

	//MpmColorTestCUDA();
	//MpmGetMomentumCUDA();

	TransferFromCUDA(fbuf);
}

void FluidSystem::RunMPM() {
	ClearTimer();

	//sorting particles
	GetParticleIndexCUDA();
	GetGridListCUDA();
	RearrageDataCUDA();
	CheckTimer("sorting");

	//node mass,density,momentum
	//MpmParticleToGrid_CUDA();
	MpmParticleToGrid_APIC_CUDA();

	//particle stress
	MpmParticleStress_CUDA();

	//node force, momentum update
	MpmNodeUpdate_CUDA();

	//particle update
	//MpmParticleUpdate_CUDA();
	MpmParticleUpdate_APIC_CUDA();

	TransferFromCUDA(fbuf);
}

void FluidSystem::RunIISPH() {
	ClearTimer();

	GetParticleIndexCUDA();
	GetGridListCUDA();
	RearrageDataCUDA();
	CheckTimer("sorting");

	ComputeDensityIISPH_CUDA();
	CheckTimer("density");

	PredictAdvection();
	CheckTimer("IISPH Predict Advection");

	PressureSolve();
	CheckTimer("IISPH Pressure Solve");
	
	//Integration();
	CheckTimer("IISPH Integration");

	TransferFromCUDA(fbuf);
}

void FluidSystem::RunSolid() {
	ClearTimer();

	//------------------Sorting-------------------
	GetParticleIndexCUDA();
	GetGridListCUDA();
	RearrageDataCUDA();
	CheckTimer("sorting");

	//density, pressure
	MfComputePressureCUDA();
	CheckTimer("pressure");
	
	//Strain, Stress
	ComputeSolidTensorX_CUDA();
	CheckTimer("Strain-Stress");

	//Force
	ComputeSolidForceX_CUDA();
	CheckTimer("Force");

	//Advance
	MfAdvanceCUDA();
	CheckTimer("advance");

	TransferFromCUDA(fbuf);
}


void FluidSystem::EmitParticles (int cat)
{
	int emitInterval = 22;
	int oldPointNum = pointNum;
    if (pointNum == maxPointNum)
        return;
	std::vector<cfloat3> plist;

	if (  frameNo % emitInterval == 0 ) {
		ParticleEmitter emitter;
		emitter.dir.Set(0,-1,0);
		emitter.position.Set(0,10,0);
		emitter.vel = 1;
		emitter.radius = 6;
		emitter.spacing = pSpacing;
		emitter.EmitImmediate(plist);

		emitter.position.Set(0, 10+pSpacing, 0);
		emitter.EmitImmediate(plist);

		InsertParticles(emitter, plist, cat);
		EmitUpdateCUDA(oldPointNum, pointNum, fbuf);
	}
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

			cfloat3 rand3(rand(), rand(), rand());
			rand3 = rand3/RAND_MAX*spacing - spacing*0.5;

			cfloat3 tmp( i*spacing, 0, j*spacing);
			rand3.y = 0;
			tmp += rand3;

			if( dot(tmp,tmp)>radius*radius)
				continue;
			mvprod(rotationMat, tmp, tmp);
			tmp += position;
			plist.push_back(tmp);
		}
	}
}

void FluidSystem::InsertParticles(ParticleEmitter& emitter, std::vector<cfloat3> plist, int cat){
	for(int i=0; i<plist.size();i++){
		int p = AddParticle();
		if(p<0)
			break;
		
		displayBuffer[p].pos = plist[i];
		//displayBuffer[p].color.Set((float)rand()/RAND_MAX, (float)rand()/RAND_MAX, (float)rand()/RAND_MAX,1.0);
		displayBuffer[p].color.Set(frameNo/200.0, 0.5, 0.8, 0.5);
		displayBuffer[p].type = 0;

		calculationBuffer[p].vel = emitter.dir * emitter.vel;
		calculationBuffer[p].veleval = emitter.dir * emitter.vel;
		calculationBuffer[p].restdens = hostCarrier.densArr[cat];
		calculationBuffer[p].mass = hostCarrier.massArr[cat];
		calculationBuffer[p].visc = hostCarrier.viscArr[cat];
	}
}

int FluidSystem::AddParticle()
{
	if (pointNum >= maxPointNum) return -1;
	int n = pointNum;

	calculationBuffer[n].vel.Set(0, 0, 0);
	calculationBuffer[n].veleval.Set(0, 0, 0);
	calculationBuffer[n].bornid = n;
	calculationBuffer[n].deformGrad.Set(unitMatrix);
	calculationBuffer[n].pressure = 0;
	calculationBuffer[n].stress.Set(0.0f);
	calculationBuffer[n].B.Set(0.0f);

	pointNum++;
	return n;
}

void FluidSystem::AddFluidVolume(cfloat3 min, cfloat3 max, float spacing, cfloat3 offs, int cat)
{
	if (pointNum==maxPointNum)
		printf("Max pointnum reached.\n");

	cfloat3 pos;
	int n=0,p;
	int cntx, cnty, cntz;
	cfloat3 cnts = (max - min) / spacing;
	cfloat3 centre = min + (max-min)/2.0f;

	float randx[3];
	float ranlen = 0.4;

	for (int y = 0; y < cnts.y; y ++) {
		for (int z=0; z < cnts.z; z++) {
			for (int x=0; x < cnts.x; x++) {
				cfloat3 rand3(rand(), rand(), rand());
				rand3 = rand3/RAND_MAX*ranlen - ranlen*0.5;
				pos = cfloat3(x, y, z)*spacing + min + rand3;

				p = AddParticle();
				if (p >= 0) {
					displayBuffer[p].pos = pos+offs;
					//displayBuffer[p].color.Set((float)rand()/RAND_MAX, (float)rand()/RAND_MAX, (float)rand()/RAND_MAX,1.0);
					displayBuffer[p].color.Set(0.2, 0.5, 0.8, 0.5);
					displayBuffer[p].type = TYPE_FLUID;

					calculationBuffer[p].restdens = hostCarrier.densArr[cat];
					calculationBuffer[p].mass = hostCarrier.massArr[cat];
					calculationBuffer[p].visc = hostCarrier.viscArr[cat];
					calculationBuffer[p].vel = cfloat3(0,-5,0);
					calculationBuffer[p].solidG = solidG;
					calculationBuffer[p].solidK = solidK;
					//calculationBuffer[p].vel.Set(centre.y-pos.y, pos.x-centre.x, 0);
					//calculationBuffer[p].vel *= 0.1;
					n++;
				}
				//return;
			}
		}
	}
	printf("No.%d fluid has %d particles\n", cat, n);
}

void FluidSystem::AddGranularVolume(cfloat3 min, cfloat3 max, float spacing, cfloat3 offs)
{
	if (pointNum==maxPointNum)
		printf("Max pointnum reached.\n");

	cfloat3 pos;
	int n=0, p;
	int cntx, cnty, cntz;
	cfloat3 cnts = (max - min) / spacing;
	cfloat3 centre = min + (max-min)/2.0f;

	float randx[3];
	float ranlen = 0.4;

	for (int y = 0; y < cnts.y; y ++) {
		for (int z=0; z < cnts.z; z++) {
			for (int x=0; x < cnts.x; x++) {
				cfloat3 rand3(rand(), rand(), rand());
				rand3 = rand3/RAND_MAX*ranlen - ranlen*0.5;
				pos = cfloat3(x, y, z)*spacing + min + rand3;

				p = AddParticle();
				if (p >= 0) {
					displayBuffer[p].pos = pos+offs;
					displayBuffer[p].color.Set(0.7+(float)rand()/4/RAND_MAX, 0.6+(float)rand()/4/RAND_MAX, 0.1+(float)rand()/4/RAND_MAX,1.0);
					//displayBuffer[p].color.Set(0.2, 0.5, 0.8, 0.5);
					displayBuffer[p].type = TYPE_GRANULAR;

					calculationBuffer[p].restdens = hostCarrier.densArr[2];
					calculationBuffer[p].mass = hostCarrier.massArr[2];
					calculationBuffer[p].visc = hostCarrier.viscArr[2];
					calculationBuffer[p].solidG = solidG;
					calculationBuffer[p].solidK = solidK;
					//calculationBuffer[p].vel = cfloat3(1,2,3);
					//calculationBuffer[p].vel.Set(centre.y-pos.y, pos.x-centre.x, 0);
					//calculationBuffer[p].vel *= 0.1;
					n++;
				}
				//return;
			}
		}
	}
	printf(" granular has %d particles\n", n);
}

void FluidSystem::AddShape(BI2Reader bi2reader, int cat)
{
	float x, y, z;

	for (int i = 0; i<bi2reader.info.np; i++)
	{
		int p = AddParticle();

		if (p!=-1)
		{
			cfloat3 offs(0,0,0);
			float scale = 0.005;
			//	printf("%f %f %f\n",bi2reader.info.Pos[i].x/m_Param[PSIMSCALE],bi2reader.info.Pos[i].y/m_Param[PSIMSCALE],bi2reader.info.Pos[i].z/m_Param[PSIMSCALE]);
			//(mPos+p)->Set (bi2reader.info.Pos[i].x/m_Param[PSIMSCALE],bi2reader.info.Pos[i].y/m_Param[PSIMSCALE],bi2reader.info.Pos[i].z/m_Param[PSIMSCALE]);
			x=bi2reader.info.Pos[i].x/scale;
			z=bi2reader.info.Pos[i].y/scale;
			y=bi2reader.info.Pos[i].z/scale;

			calculationBuffer[p].vel = cfloat3(0, -1.5, 0);
			displayBuffer[p].pos = cfloat3(x,y,z) + offs;
			//displayBuffer[p].color.Set(0.7+(float)rand()/4/RAND_MAX, 0.6+(float)rand()/4/RAND_MAX, 0.1+(float)rand()/4/RAND_MAX, 1.0);
			displayBuffer[p].color.Set(0.2, 0.5, 0.8, 0.5);
			displayBuffer[p].type = TYPE_FLUID;

			calculationBuffer[p].restdens = hostCarrier.densArr[2];
			calculationBuffer[p].mass = hostCarrier.massArr[2];
			calculationBuffer[p].visc = hostCarrier.viscArr[2];
			calculationBuffer[p].solidG = solidG;
			calculationBuffer[p].solidK = solidK;
		}
	}
}

void FluidSystem::AddDeformableVolume(cfloat3 min, cfloat3 max, float spacing, int cat)
{
	if (pointNum==maxPointNum)
		printf("Max pointnum reached.\n");

	cfloat3 pos;
	int n = 0, p;
	int cntx, cnty, cntz;
	cfloat3 cnts = (max - min) / spacing;
	cfloat3 center = min + (max-min)*0.5;

	float randx[3];
	float ranlen = 0.2;

	for (int y = 0; y < cnts.y; y ++) {
		for (int z=0; z < cnts.z; z++) {
			for (int x=0; x < cnts.x; x++) {
				cfloat3 rand3(rand(), rand(), rand());
				rand3 = rand3/RAND_MAX*ranlen - ranlen*0.5;
				pos = cfloat3(x, y, z)*spacing + min + rand3;

				p = AddParticle();
				if (p >= 0) {
					displayBuffer[p].pos = pos;
					displayBuffer[p].color.Set(0.8, 0.5, 0.2, 0.5);
					displayBuffer[p].type = TYPE_ELASTIC;

					calculationBuffer[p].restdens = hostCarrier.densArr[cat];
					calculationBuffer[p].mass = hostCarrier.massArr[cat];
					calculationBuffer[p].visc = hostCarrier.viscArr[cat];
					
					calculationBuffer[p].vel.y = (pos.x-center.x)/(max.x-min.x)*2.5;
					calculationBuffer[p].X = pos;
				}
			}
		}
	}
	printf("%d fluid has %d particles\n", cat, n);
}

void FluidSystem::SetDruckerPragerParam() {
	float phi = hostCarrier.a_phi;
	float psi = hostCarrier.a_psi;
	float k_c = hostCarrier.k_c;

	phi *= 3.141593f/180.0f;
	psi *= 3.141593f/180.0f;
	hostCarrier.a_phi = 2*sin(phi)/sqrt(3)/(3-sin(phi));
	hostCarrier.a_psi = 2*sin(psi)/sqrt(3)/(3-sin(psi));
	hostCarrier.k_c = 3*k_c*cos(psi)/sqrt(3)/(3-sin(phi));
}

using namespace tinyxml2;
extern XMLElement* pele;

void FluidSystem::ParseXML(int caseid){
	
	printf("Parsing XML.\n");
	tinyxml2::XMLDocument doc;
	int result = doc.LoadFile("scene.xml");
	

	XMLElement* fluidElement = doc.FirstChildElement("Fluid");
	XMLElement* sceneElement = doc.FirstChildElement("MultiScene");
	XMLElement* boundElement = doc.FirstChildElement("BoundInfo");

	if(!fluidElement || !sceneElement)
	{
		printf("missing fluid/scene xml node");
		return;
	}

	int tmp;
	
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
	hostCarrier.acclimit = XMLGetFloat("AccelLimit");
	hostCarrier.vlimit = XMLGetFloat("VelLimit");
	hostCarrier.gravity = XMLGetFloat3("Gravity");
	hostCarrier.volmin = XMLGetFloat3("VolMin");
	hostCarrier.volmax = XMLGetFloat3("VolMax");


	//scene specified
	pele = sceneElement;
	
	maxPointNum = XMLGetInt("Num");
	scaleP = XMLGetFloat("ScaleP");
	fluidnum = XMLGetInt("FluidCount");
	fvmin[0] = XMLGetFloat3("VolMin0");
	fvmax[0] = XMLGetFloat3("VolMax0");
	fvmin[1] = XMLGetFloat3("VolMin1");
	fvmax[1] = XMLGetFloat3("VolMax1");
	XMLGetFloatN(massratio, 3, "MassRatio");
	XMLGetFloatN(densratio,3,"DensityRatio");
	XMLGetFloatN(viscratio, 3, "ViscRatio");

	hostCarrier.maxNum = maxPointNum;
	hostCarrier.simscale = XMLGetFloat("SimScale");
	hostCarrier.dt = XMLGetFloat("DT");
	hostCarrier.smoothradius = XMLGetFloat("SmoothRadius");
	hostCarrier.viscosity = XMLGetFloat("Viscosity");
	hostCarrier.mass = XMLGetFloat("Mass");
	hostCarrier.restdensity = XMLGetFloat("RestDensity");
	hostCarrier.radius = XMLGetFloat("Radius");
	hostCarrier.softminx = XMLGetFloat3("SoftBoundMin");
	hostCarrier.softmaxx = XMLGetFloat3("SoftBoundMax");
	hostCarrier.solidK = XMLGetFloat("SolidK");
	hostCarrier.solidG = XMLGetFloat("SolidG");
	hostCarrier.boundstiff = XMLGetFloat("BoundStiff");
	
	hostCarrier.a_phi = XMLGetFloat("DP_Phi");
	hostCarrier.a_psi = XMLGetFloat("DP_Psi");
	hostCarrier.k_c = XMLGetFloat("DP_cohesion");
	SetDruckerPragerParam();

	//boundary
	pele = boundElement;
	hostCarrier.bmass = XMLGetFloat("Mass");
	hostCarrier.bRestdensity = XMLGetFloat("RestDensity");
	hostCarrier.bvisc = XMLGetFloat("Viscosity");

}

void FluidSystem::TransferToCUDA(bufList& fbuf) {
	cudaMemcpy(fbuf.displayBuffer,displayBuffer,sizeof(displayPack)*pointNum, cudaMemcpyHostToDevice);
	cudaMemcpy(fbuf.calcBuffer, calculationBuffer, sizeof(calculationPack)*pointNum, cudaMemcpyHostToDevice);
}

void FluidSystem::TransferFromCUDA(bufList& fbuf) {
	cudaMemcpy( displayBuffer, fbuf.displayBuffer, sizeof(displayPack)*pointNum, cudaMemcpyDeviceToHost);
	//cudaMemcpy( calculationBuffer, fbuf.calcBuffer, sizeof(calculationPack)*pointNum, cudaMemcpyDeviceToHost);
}


void FluidSystem::EmitUpdateCUDA(int startnum, int endnum, bufList& fbuf)
{
	cudaMemcpy(fbuf.displayBuffer+startnum, displayBuffer+startnum, sizeof(displayPack)*(endnum-startnum), cudaMemcpyHostToDevice);
	cudaMemcpy(fbuf.calcBuffer+startnum, calculationBuffer+startnum, sizeof(calculationPack)*(endnum-startnum), cudaMemcpyHostToDevice);

	hostCarrier.num = pointNum;
	CarryParam(hostCarrier);
	fcuda.pnum = pointNum;
	computeNumBlocks(fcuda.pnum, 384, fcuda.numBlocks, fcuda.numThreads);
}


//-----------------	MPM Grid -----------------

void FluidSystem::SetupMPMGrid() {
	//parameter
	mpmxmin = hostCarrier.softminx - cfloat3(10,10,10);
	mpmxmax = hostCarrier.softmaxx + cfloat3(10,10,10);
	mpmcellsize = pSpacing * 1.5;
	cfloat3 tmp = (mpmxmax - mpmxmin) / mpmcellsize + cfloat3(2,2,2); //node number of each dimension
	mpmres.Set(tmp.x,tmp.y,tmp.z);
	mpmxmax = mpmxmin + (cfloat3)(mpmres+cint3(-1,-1,-1)) * mpmcellsize;

	hostCarrier.mpmXmin = mpmxmin;
	hostCarrier.mpmXmax = mpmxmax;
	hostCarrier.mpmcellsize = mpmcellsize;
	hostCarrier.mpmNodeNum = mpmres.prod();
	hostCarrier.mpmRes = mpmres;

	nodevel = new cfloat3[ hostCarrier.mpmNodeNum ];

	MpmAllocateBufferCUDA(hostCarrier);
	CarryParam(hostCarrier);

	IndexMPMSortCUDA();
}

void FluidSystem::ReleaseMPMGrid() {
	delete nodevel;
}

void FluidSystem::CopyMPMFromDevice() {
	cudaMemcpy(nodevel, fbuf.mpmVel, sizeof(cfloat3)*hostCarrier.mpmNodeNum, cudaMemcpyDeviceToHost);
}

//--------------------End MPM Grid-----------------


//--------------------Solid------------------------

