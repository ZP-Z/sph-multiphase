

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

#include <conio.h>

//#include <GL/glut.h>
//#include <cuda_gl_interop.h>

#include "fluid_system_host.cuh"		
#include "fluid_system_kern.cuh"
//#include "radixsort.cu"						// Build in RadixSort
#include "thrust\device_vector.h"	//thrust libs
#include "thrust\sort.h"
#include "fluidMath.cuh"

__device__ ParamCarrier paramCarrier;


void CarryParam(ParamCarrier& hostCarrier){
	cudaMemcpyToSymbol( paramCarrier, &hostCarrier, sizeof(hostCarrier));
}

void updateParam(FluidParams* paramCPU){
    //cudaMemcpyToSymbol ( simData, paramCPU, sizeof(FluidParams) );
}



//Sorting

//get particle index id
__global__ void InitialSort ( bufList buf, int pnum )
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;
	register cfloat3 gridMin = paramCarrier.gridmin;
	register cfloat3 gridDelta = paramCarrier.gridIdfac;
	register cint3 gridRes = paramCarrier.gridres;
	register cint3 gridScan = paramCarrier.gridres;
	gridScan.x -= 1;
	gridScan.y -= 1;
	gridScan.z -= 1;

	register int		gs;
	register cfloat3		gcf;
	register cint3		gc;

	gcf = (buf.displayBuffer[i].pos - gridMin) * gridDelta; 
	gc = cint3( int(gcf.x), int(gcf.y), int(gcf.z) );
	
	gs = (gc.y * gridRes.z + gc.z)*gridRes.x + gc.x;
	if ( gc.x >= 1 && gc.x < gridScan.x && gc.y >= 1 && gc.y < gridScan.y && gc.z >= 1 && gc.z < gridScan.z ) {
		buf.mgcell[i] = gs;											// Grid cell insert.
		buf.midsort[i] = i;
	} else {
		buf.mgcell[i] = GRID_UNDEF;
		buf.midsort[i] = i;
	}

}

//markup the head and tail of each cell
__global__ void CalcFirstCnt ( bufList buf, int pnum )
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i>=pnum) return;
	if ((i == 0 || buf.mgcell[i]!=buf.mgcell[i-1]))
	{
		if (buf.mgcell[i]!=GRID_UNDEF)
			buf.mgridoff[buf.mgcell[i]] = i; //head id - gridoff
	}
	
	if (i!=0 && buf.mgcell[i]!=buf.mgcell[i-1] && buf.mgcell[i-1]!=GRID_UNDEF)
		buf.mgridcnt[buf.mgcell[i-1]] = i; //tail id + 1 - gridcnt
	
	if (i == pnum-1 && buf.mgcell[i]!=GRID_UNDEF)
		buf.mgridcnt[buf.mgcell[i]] = i + 1;	
}

__global__ void GetCnt ( bufList buf, int pnum )
{
	uint i = blockIdx.x * blockDim.x + threadIdx.x;		// particle index
	if (i>=pnum) return ;
	if (buf.mgcell[i]!=GRID_UNDEF)
	{
		buf.mgndx[i] = i - buf.mgridoff[buf.mgcell[i]];

		if (buf.mgndx[i] == 0) // first particle of the grid, once
			buf.mgridcnt[buf.mgcell[i]] -= buf.mgridoff[buf.mgcell[i]]; //cnt = tail - head
	}
}

//deep copy sort
__global__ void RearrangeData( bufList buf, int pnum )
{
	//for each new position, find old particle and read value
	uint i = blockIdx.x * blockDim.x + threadIdx.x;		// particle index				
	if ( i >= pnum ) return;

	int exId = buf.midsort[i]; //original id
	int cell = buf.mgcell[i];
	int sort_ndx = i;
	i = exId;

	
	if ( cell != GRID_UNDEF ) {
		
		buf.MFidTable[ buf.calcBuffer[sort_ndx].bornid ] = sort_ndx;
		
		buf.displayBuffer[sort_ndx] = *(displayPack*)(buf.msortbuf+pnum*BUF_DISPLAYBUF+i*sizeof(displayPack));
		buf.calcBuffer[sort_ndx] = *(calculationPack*)(buf.msortbuf+pnum*BUF_CALCBUF+i*sizeof(calculationPack));
	}
	else{
		buf.mgcell[sort_ndx] = GRID_UNDEF;
		buf.displayBuffer[sort_ndx].pos.Set(-1000,-1000,-1000);
	}
}

//__device__ void findNearest(int i, float& mindis, int cell, bufList buf)
//{
//	cfloat3 dist;
//	float dsq;
//	cfloat3 p = buf.displayBuffer[i].pos;
//
//	register float d2 = paramCarrier.simscale * paramCarrier.simscale;
//	register float r2 = simData.r2/d2;
//	//int j;
//
//	if (buf.mgridcnt[cell] == 0) return;
//
//	int cfirst = buf.mgridoff[cell];
//	int clast = cfirst + buf.mgridcnt[cell];
//	for (int j = cfirst; j < clast; j++) {
//		//j = buf.mgrid[cndx];
//		if (buf.displayBuffer[i].type != TYPE_BOUNDARY) {
//			dist = p - buf.displayBuffer[j].pos;
//			dsq = dot(dist, dist);
//
//			if (dsq < r2 && dsq > 0.0 && dsq*d2<mindis)
//			{
//				mindis = dsq*d2;
//				buf.midsort[i] = j;
//			}
//		}
//	}
//
//	return;
//}
//__global__ void mfFindNearest (bufList buf,int pnum)
//{
//	uint i = blockIdx.x * blockDim.x + threadIdx.x;				
//	if ( i >= pnum ) return;
//	
//	// Get search cell
//	uint gc = buf.mgcell[i];
//	if ( gc == GRID_UNDEF ) return;
//	
//	// Sum Pressures
//	cfloat3 pos = buf.displayBuffer[i].pos;
//	float mindis = 65535;
//
//	if (buf.displayBuffer[i].type == TYPE_BOUNDARY) //boundary particles
//	{
//		buf.midsort[i] = i;
//		buf.calcBuffer[i].mass = simData.pmass;
//		for (int c = 0; c<simData.gridAdjCnt; c++)
//		{
//			findNearest(i,mindis,gc+simData.gridAdj[c],buf);
//		}
//		
//		///nearest id ---> idsort
//
//		if (buf.midsort[i]!=i)
//			buf.calcBuffer[i].mass  = buf.calcBuffer[buf.midsort[i]].mass;
//
//	}
//}

__device__ void contributeDensity_boundary(uint i, float& res, uint cell, bufList buf) {
	cfloat3 dist;
	cfloat3 p = buf.displayBuffer[i].pos;
	float dsq, c;
	float massj;
	float r2 = paramCarrier.smoothradius;
	r2 = r2*r2;

	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	for (int j = cfirst; j < clast; j++) {

		if ( buf.displayBuffer[j].type != TYPE_BOUNDARY )
			continue;

		dist = p - buf.displayBuffer[j].pos;
		dist = dist * paramCarrier.simscale;

		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if (dsq < r2) {
			c = r2 - dsq;
			res += c * c * c ;//* buf.mf_restmass[j];
		}
	}
	return;
}

__device__ void contributeDensity (uint i, float& res, uint cell, bufList buf){
	cfloat3 dist;
	cfloat3 p = buf.displayBuffer[i].pos;
	float dsq, c;
	float massj;
	float r2 = paramCarrier.smoothradius;
	r2 = r2*r2;

	//int j;

	if ( buf.mgridcnt[cell] == 0 )
		return;

	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];

	for ( int j = cfirst; j < clast; j++ ) {
		
        if ( buf.displayBuffer[j].type==0) //fluid
			massj = buf.calcBuffer[j].mass;
		else if(buf.displayBuffer[j].type==1) //boundary
			massj = buf.calcBuffer[i].restdens * buf.calcBuffer[j].dens; // fluid density * boundary volume

		dist = p - buf.displayBuffer[j].pos;
		dist = dist * paramCarrier.simscale;

		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if ( dsq < r2 ) {
			c = r2 - dsq;
			res += c * c * c * massj;	
		} 
	}
	return;
}

__global__ void ComputeBoundaryVolume(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if (i >= pnum) return;

	if(buf.displayBuffer[i].type != 1) //not boundary
		return;

	// Get search cell
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	

	float sum = 0.0;

	//Get Boundary Density
	
	for (int c=0; c < paramCarrier.neighbornum; c++) {
		contributeDensity_boundary(i, sum, gc + paramCarrier.neighborid[c], buf);
	}

	sum = sum * paramCarrier.kpoly6;

	if (sum == 0.0) {
		printf("boundary density zero error.\n");
		sum = 1.0;
	}
		
	buf.calcBuffer[i].dens = 1 / sum ; //actually the volume
}

__global__ void ComputeDensityPressure(bufList buf,int pnum){
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	if(buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;

	// Get search cell
	uint gc = buf.mgcell[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	

	float sum = 0.0;
	float dens;
	
	dens = buf.calcBuffer[i].restdens;

	//Get Fluid Density
	for (int c=0; c < paramCarrier.neighbornum; c++) {
		contributeDensity(i, sum, gc + paramCarrier.neighborid[c], buf);
	}

	sum = sum * paramCarrier.kpoly6;
	
	if ( sum == 0.0 )
		sum = 1.0;
	
	buf.calcBuffer[i].dens = 1/sum;

	//buf.mpress[i] = ( sum - dens ) * simData.pintstiff;
	//buf.mpress[i] = (pow( sum/dens,7.0f )-1) * simData.pintstiff;
    buf.calcBuffer[i].pressure = 2.5 * dens * (pow(sum / dens, 7.0f) - 1);

    if (buf.calcBuffer[i].pressure<0)
        buf.calcBuffer[i].pressure = 0;
	//if(buf.calcBuffer[i].bornid == 0)
	//	printf("%f\n",buf.calcBuffer[i].pressure);

	//buf.mclr[i] = COLORA(1, 1-buf.mpress[i]/1000, 1-buf.mpress[i]/1000, 1);
}






__device__ cfloat3 contributeForce_new(int i, cfloat3 ipos, cfloat3 iveleval, float ipress, float idens, int cell, bufList buf, cfloat3* ivelxcor, float ivisc)
{
	//Force here represents the acceleration
	float dsq, c;
	float r2 = paramCarrier.smoothradius;
	r2 = r2*r2;

	cfloat3 dist, vmr;
	float cmterm, cmterm1;

	float pmterm, vmterm;

	int j;
	float aveDenij, cx, xterm;

	if (buf.mgridcnt[cell] == 0) return cfloat3(0, 0, 0);

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	cfloat3 force = cfloat3(0, 0, 0);

	for (int j = cfirst; j < clast; j++) {
		

		dist = (ipos - buf.displayBuffer[j].pos);
		dist *= paramCarrier.simscale;
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		
		if (!(dsq < r2 && dsq > 0))
			continue;

		cx = (r2-dsq);
		dsq = sqrt(dsq);
		c = (paramCarrier.smoothradius - dsq);

		cmterm1 = paramCarrier.kspikydiff * c * c / dsq; //nabla W
		cmterm = cmterm1 * buf.calcBuffer[j].mass * buf.calcBuffer[j].dens;


		if (buf.displayBuffer[j].type == 0) 
		{
			//pressure
			pmterm = - 0.5* cmterm * (buf.calcBuffer[j].pressure) *idens;
			force += dist * pmterm;

			//viscosity
			vmr = iveleval - buf.calcBuffer[j].vel;
			vmterm = cmterm * (ivisc+buf.calcBuffer[j].visc) * idens;
			force += vmr*vmterm;

		}
		else if(buf.displayBuffer[j].type == 1){
			//pressure
			pmterm = - cmterm1 * buf.calcBuffer[i].restdens * buf.calcBuffer[j].dens *  (ipress ) *idens*idens;
			force += dist * pmterm;
			
			//artificial boundary viscosity			
			vmr = iveleval - buf.calcBuffer[j].vel;
			float pi_ij = vmr.x*dist.x + vmr.y*dist.y + vmr.z*dist.z;
			if (pi_ij < 0) {
				pi_ij = pi_ij / (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z + r2 * 0.01);
				pi_ij = pi_ij * paramCarrier.bvisc * paramCarrier.smoothradius  /2 / buf.calcBuffer[i].restdens ;
				pi_ij = cmterm1 * buf.calcBuffer[i].restdens * buf.calcBuffer[j].dens * pi_ij;
				force += (dist * pi_ij);
			}
		}
		
	}
	return force;
}

__global__ void ComputeForce (bufList buf, int pnum){
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum)
		return;
	if(buf.displayBuffer[i].type==1)
		return;

	// Get search cell
	uint gc = buf.mgcell[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	

	register cfloat3 ipos = buf.displayBuffer[i].pos;
	register cfloat3 iveleval = buf.calcBuffer[i].vel;
	register float ipress = buf.calcBuffer[i].pressure;
	register float idens = buf.calcBuffer[i].dens;
	//register cfloat3 ivelxcor = buf.calcBuffer[i].velxcor;
	register float ivisc = buf.calcBuffer[i].visc;

	register cfloat3 force = cfloat3(0,0,0);	
	cfloat3 ivelxcor = cfloat3(0,0,0);

	for (int c=0; c < paramCarrier.neighbornum; c++) {
		force += contributeForce_new ( i,  ipos, iveleval, ipress, idens, gc + paramCarrier.neighborid[c], buf, &ivelxcor, ivisc);
	}
	buf.calcBuffer[i].accel = force;

}

__global__ void AdvanceParticles(bufList buf, int numPnts)
{

	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if (i >= numPnts) return;
	if (buf.displayBuffer[i].type==1)
		return;

	if (buf.mgcell[i] == GRID_UNDEF) {
		buf.displayBuffer[i].pos = cfloat3(-1000, -1000, -1000);
		buf.calcBuffer[i].vel = cfloat3(0, 0, 0);
		return;
	}

	// Get particle vars
	cfloat3 accel, norm;
	float diff, adj, accelLen;
	cfloat3 pos = buf.displayBuffer[i].pos;
	cfloat3 veval = buf.calcBuffer[i].vel;
	float dt = paramCarrier.dt;
	float ss = paramCarrier.simscale;

	accel = cfloat3(0, 0, 0);
	
	/*
	// Soft Boundaries
	// Y-axis
	diff = simData.pradius - (pos.y - (simData.pboundmin.y + (pos.x-simData.pboundmin.x)*simData.pground_slope)) * ss;
	float slope = 0;
	if (diff > EPSILON) {
		norm = cfloat3(-slope, 1.0 -slope, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;
		//norm = -cfloat3(vel.x, 0, vel.z);
		//accel += norm*10;

	}

	diff = simData.pradius - (simData.pboundmax.y - pos.y)*ss;
	if (diff > EPSILON) {
		norm = cfloat3(0, -1, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;
	}

	// X-axis
	diff = simData.pradius - (pos.x - (simData.pboundmin.x + (sin(time*simData.pforce_freq)+1)*0.5 * simData.pforce_min))*ss;
	if (diff > EPSILON) {
		norm = cfloat3(1, 0, 0);
		adj = (simData.pforce_min+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}

	diff = simData.pradius - ((simData.pboundmax.x - (sin(time*simData.pforce_freq)+1)*0.5*simData.pforce_max) - pos.x)*ss;
	if (diff > EPSILON) {
		norm = cfloat3(-1, 0, 0);
		adj = (simData.pforce_max+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}

	// Z-axis
	diff = simData.pradius - (pos.z - simData.pboundmin.z) * ss;
	if (diff > EPSILON) {
		norm = cfloat3(0, 0, 1);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}

	diff = simData.pradius - (simData.pboundmax.z - pos.z)*ss;
	if (diff > EPSILON) {
		norm = cfloat3(0, 0, -1);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}
	*/
	//End Soft Boundary

	accel += buf.calcBuffer[i].accel;
	accel += paramCarrier.gravity;

	// Accel Limit
	//accelLen = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
	//if (accelLen > simData.AL2) {
	//	accel *= simData.AL / sqrt(accelLen);
	//}

	cfloat3 vel = buf.calcBuffer[i].vel;
	//Velocity Limit
	/*float vmod = sqrtf(dot(vel, vel));
	if (vmod > simData.fluidVConstraint)
	vel *= simData.fluidVConstraint / vmod;*/

	cfloat3 vnext = accel*dt + vel;		// v(t+1/2) = v(t-1/2) + a(t) dt
	buf.displayBuffer[i].pos += vnext*dt/ss;
	buf.calcBuffer[i].veleval = (vel+vnext)*0.5;
	buf.calcBuffer[i].vel = vnext;

	//For Drift Velocity Calculation of Next Step
	//buf.mforce[i] = simData.pgravity - accel;
	//buf.mforce[i] = accel;
	//buf.mforce[i] = cfloat3(0, 0, 0);

	//Color Setting
	/*if(buf.MFtype[i]==3)
	buf.mclr[i] = COLORA(buf.mf_alpha[i*MAX_FLUIDNUM+2],buf.mf_alpha[i*MAX_FLUIDNUM+1],buf.mf_alpha[i*MAX_FLUIDNUM+0],0.6);
	else
	buf.mclr[i] = COLORA(buf.mf_alpha[i*MAX_FLUIDNUM+2],buf.mf_alpha[i*MAX_FLUIDNUM+1],buf.mf_alpha[i*MAX_FLUIDNUM+0],1);*/

}





/*-------------------------

	Surface Tension - Yang Tao 

--------------------------*/



__device__ void contributeCover(int i, int cell, bufList buf){

}

__global__ void SurfaceDetection(bufList buf, int pnum){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum)
		return;
	
	// Get search cell
	uint gc = buf.mgcell[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	

	for (int c=0; c < paramCarrier.neighbornum; c++) {
		contributeCover ( i, gc + paramCarrier.neighborid[c], buf );
	}
}



__device__ void contributeSurfaceTensionYT(int i, cfloat3& res, int cell, bufList buf)
{
    //Force here represents the acceleration
	float dsq, c;
	float sr2 = paramCarrier.smoothradius;
	sr2 *= sr2;
	cfloat3 dist;		
	int j;	

	if ( buf.mgridcnt[cell] == 0 )
        return;	

	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];

    float tensionRadius2 = paramCarrier.surfaceTensionK * paramCarrier.surfaceTensionK * sr2;
    float tR = sqrt(tensionRadius2);

    float forceMod;
    float ffTensionC = paramCarrier.surfaceTensionFluidC;

	for ( int j = cfirst; j < clast; j++ ) {	

        dist = buf.displayBuffer[i].pos - buf.displayBuffer[j].pos;
        dist = dist * paramCarrier.simscale;
		dsq = dot(dist,dist);
        if (dsq != 0.0 && dsq <= tR){
            dsq = sqrt(dsq);
            forceMod = ffTensionC * cos( 1.5 * 3.14159 * dsq/ tR) / dsq;
            
            res += dist*forceMod;
        }

    }
}

__global__ void SurfaceTension(bufList buf, int pnum){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum)
		return;
	//if(buf.displayBuffer[i].type==1)
	//	return;

	// Get search cell
	uint gc = buf.mgcell[i];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	

	register cfloat3 force = cfloat3(0,0,0);	

	for (int c=0; c < paramCarrier.neighbornum; c++) {
		contributeSurfaceTensionYT ( i, force, gc + paramCarrier.neighborid[c], buf);
	}

	buf.calcBuffer[i].accel += force;

}


/*-------------------------

	End Surface Tension

--------------------------*/

__device__ void contributeAcc_adv(int i, cfloat3& acc_adv, int cell, bufList buf) {
	
	if (buf.mgridcnt[cell] == 0)
		return;
	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	float c;
	float h = paramCarrier.smoothradius;
	float h2 = h*h;

	cfloat3 xij, vij;
	float dist;
	int j;

	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 ivel = buf.calcBuffer[i].vel;
	float ivisc = buf.calcBuffer[i].visc;
	float idens = buf.calcBuffer[i].dens;
	float cmterm;
	float nW_fac;
	float vmterm;
	cfloat3 accel = cfloat3(0, 0, 0);

	for (int j = cfirst; j < clast; j++) {

		xij = ipos - buf.displayBuffer[j].pos;
		xij *= paramCarrier.simscale;
		dist = sqrt(dot(xij,xij));

		if (!(dist < h && dist > 0))
			continue;
		c = h - dist;

		nW_fac = paramCarrier.kspikydiff * c * c / dist;
		cmterm = nW_fac * buf.calcBuffer[j].mass * buf.calcBuffer[j].dens;

		if (buf.displayBuffer[j].type == 0)
		{
			vij = ivel - buf.calcBuffer[j].vel;
			vmterm = cmterm * (ivisc+buf.calcBuffer[j].visc) * idens;
			accel += vij*vmterm;
		}
		else if (buf.displayBuffer[j].type == 1) {

			//artificial boundary viscosity			
			vij = ivel - buf.calcBuffer[j].vel;
			float pi_ij = dot(xij,vij);
			if (pi_ij < 0) {
				pi_ij = pi_ij / (dist*dist + h2 * 0.01);
				pi_ij = pi_ij * paramCarrier.bvisc * paramCarrier.smoothradius  /2 / buf.calcBuffer[i].restdens;
				pi_ij = nW_fac * buf.calcBuffer[i].restdens * buf.calcBuffer[j].dens * pi_ij;
				accel += (xij * pi_ij);
			}
		}

	}
	
	acc_adv = acc_adv + accel;
}

__device__ void contributeDii(int i, cfloat3& dii, int cell, bufList buf) {
	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 jpos,xij;
	float dist;
	float h = paramCarrier.smoothradius;
	float c;
	float nW_fac;
	cfloat3 nablaWij;

	for (int j = cfirst; j < clast; j++) {
		
		jpos = buf.displayBuffer[j].pos;
		xij = ipos - jpos;
		xij = xij*paramCarrier.simscale;
		dist = sqrt(dot(xij,xij));
		if(!(dist<h && dist>0))
			continue;
		
		c = h - dist;
		nW_fac = paramCarrier.kspikydiff * c * c / dist;
		nablaWij = xij * nW_fac;
		
		if (buf.displayBuffer[j].type!=TYPE_BOUNDARY)
			dii = dii - nablaWij*buf.calcBuffer[j].mass;
		else
			dii = dii - nablaWij*buf.calcBuffer[i].restdens*buf.calcBuffer[j].dens;
	}
}

__global__ void ComputeDii(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;			
	if (i >= pnum) return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;



	//v_adv
	cfloat3 acc_adv(0,0,0);
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		contributeAcc_adv(i, acc_adv, gc+paramCarrier.neighborid[c],buf);
	}
	acc_adv += paramCarrier.gravity;
	buf.calcBuffer[i].veleval = buf.calcBuffer[i].vel + acc_adv * paramCarrier.dt;



	//dii
	cfloat3 dii(0,0,0);
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		contributeDii(i,dii,gc+paramCarrier.neighborid[c],buf);
	}
	
	dii = dii * buf.calcBuffer[i].dens * buf.calcBuffer[i].dens * paramCarrier.dt * paramCarrier.dt;
	//if (i%100==0)
	//	printf("dii: %.20f %.20f %.20f\n", dii.x, dii.y, dii.z);

	buf.dii[i] = dii;
	
}


__device__ void contributeDensChange(int i, float& denschange, int cell, bufList buf) {
	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 ivel = buf.calcBuffer[i].veleval;
	cfloat3 jpos, xij, vij;
	float dist;
	float h = paramCarrier.smoothradius;
	float c;
	float nW_fac;
	cfloat3 nablaWij;
	float sum = 0;

	for (int j = cfirst; j < clast; j++) {
		jpos = buf.displayBuffer[j].pos;
		xij = ipos - jpos;
		xij = xij*paramCarrier.simscale;
		dist = sqrt(dot(xij,xij));
		if (!(dist<h && dist>0))
			continue;

		vij = ivel - buf.calcBuffer[j].veleval;
		c = h - dist;
		nW_fac = paramCarrier.kspikydiff * c * c / dist; //nabla W
		nablaWij = xij*nW_fac;

		if(buf.displayBuffer[j].type!=TYPE_BOUNDARY)
			sum += dot(vij, nablaWij) * buf.calcBuffer[j].mass * paramCarrier.dt;
		else
			sum += dot(vij, nablaWij) * buf.calcBuffer[j].dens * buf.calcBuffer[i].restdens * paramCarrier.dt;
	}
	denschange += sum;
}

__device__ void contributeAii(int i, float& aii, int cell, bufList buf) {
	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 ivel = buf.calcBuffer[i].veleval;
	cfloat3 jpos, xij;
	float dist;
	float h = paramCarrier.smoothradius;
	float c;
	float nW_fac;
	cfloat3 nablaWij;

	float sum = 0;
	cfloat3 dji;
	cfloat3 dii = buf.dii[i];
	cfloat3 dtmp;

	for (int j = cfirst; j < clast; j++) {
		if (buf.displayBuffer[j].type==TYPE_BOUNDARY)
			continue;
		jpos = buf.displayBuffer[j].pos;
		xij = ipos - jpos;
		xij = xij*paramCarrier.simscale;
		dist = sqrt(dot(xij,xij));

		if (!(dist<h && dist>0))
			continue;
		c = h - dist;
		nW_fac = paramCarrier.kspikydiff * c * c / dist; //nabla W
		nablaWij = xij*nW_fac;

		dji = nablaWij * buf.calcBuffer[i].mass* buf.calcBuffer[i].dens * buf.calcBuffer[i].dens * paramCarrier.dt
			* paramCarrier.dt;
		dtmp = dii - dji;
		sum += dot(dtmp, nablaWij) * buf.calcBuffer[j].mass;
	}
	aii += sum;

}

__global__ void ComputeAii(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum) return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;

	//rho_i_adv
	float denschange = 0;
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		contributeDensChange(i, denschange, gc+paramCarrier.neighborid[c],buf);
	}
	buf.rho_adv[i] = 1/buf.calcBuffer[i].dens + denschange;

	//p_i_l0
	buf.press_l[i] = 0.5 * buf.calcBuffer[i].pressure;

	//aii
	float aii = 0;
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		contributeAii(i, aii, gc+paramCarrier.neighborid[c], buf);
	}
	buf.aii[i] = aii;

}

__device__ void contributeDP(int i, cfloat3& dp, int cell, bufList buf) {
	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 jpos, xij;
	float dist;
	float h = paramCarrier.smoothradius;
	float c;
	float nW_fac;
	cfloat3 nablaWij;

	for (int j = cfirst; j < clast; j++) {
		if (buf.displayBuffer[j].type==TYPE_BOUNDARY)
			continue;
		jpos = buf.displayBuffer[j].pos;
		xij = ipos - jpos;
		xij = xij*paramCarrier.simscale;
		dist = sqrt(xij.x*xij.x + xij.y*xij.y + xij.z*xij.z);
		if (!(dist<h && dist>0))
			continue;
		c = h - dist;
		nW_fac = paramCarrier.kspikydiff * c * c / dist; //nabla W
		nablaWij = xij*nW_fac;
		
		dp = dp + nablaWij * buf.press_l[j] * buf.calcBuffer[j].mass * buf.calcBuffer[j].dens
			* buf.calcBuffer[j].dens * (-1);
	}

}

__global__ void Pressure_DP(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum) return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;

	cfloat3 dp(0,0,0);
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		contributeDP(i,dp,gc+paramCarrier.neighborid[c],buf);
	}
	dp = dp * paramCarrier.dt * paramCarrier.dt;
	buf.dijpj[i] = dp;
}

__device__ void contributeP_updateterm(int i, float& updateterm, int cell, bufList buf) {
	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];
	cfloat3 ipos = buf.displayBuffer[i].pos;
	cfloat3 jpos, xij;
	float dist;
	float h = paramCarrier.smoothradius;
	float c;
	float nW_fac;
	cfloat3 nablaWij;
	
	cfloat3 dijpj = buf.dijpj[i];
	cfloat3 djj;
	cfloat3 djkpk, dji;
	cfloat3 tmpf3;

	for (int j = cfirst; j < clast; j++) {
		
		jpos = buf.displayBuffer[j].pos;
		xij = ipos - jpos;
		xij = xij*paramCarrier.simscale;
		dist = sqrt(xij.x*xij.x + xij.y*xij.y + xij.z*xij.z);
		if (!(dist<h && dist>0))
			continue;
		c = h - dist;
		nW_fac = paramCarrier.kspikydiff * c * c / dist; //nabla W
		nablaWij = xij*nW_fac;

		if (buf.displayBuffer[j].type!=TYPE_BOUNDARY) {
			djj = buf.dii[j];
			dji = nablaWij*buf.calcBuffer[i].mass*buf.calcBuffer[i].dens*buf.calcBuffer[i].dens*
				paramCarrier.dt * paramCarrier.dt;
			djkpk = buf.dijpj[j] - dji * buf.press_l[i];
			tmpf3 = dijpj - djj*buf.press_l[j] - djkpk;
			updateterm += buf.calcBuffer[j].mass * dot(tmpf3, nablaWij);
		}
		else {
			updateterm += buf.calcBuffer[i].restdens * buf.calcBuffer[j].dens * dot(dijpj,nablaWij);
		}
		
	}

}


__global__ void Pressure_Iter(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum) return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY){
		buf.densityResidue[i] = 0;
		return;
	}
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;

	//compute update term
	float updateterm = 0;
	for (int c=0; c<paramCarrier.neighbornum; c++) {
		contributeP_updateterm(i, updateterm, gc+paramCarrier.neighborid[c], buf);
	}

	//update pressure
	float omega = 0.3;
	if (abs(buf.aii[i])<EPSILON) {
		buf.press_l1[i] = buf.press_l[i];
	}
	else{
		buf.press_l1[i] = (1-omega)*buf.press_l[i] + omega/buf.aii[i]
			*(buf.calcBuffer[i].restdens - buf.rho_adv[i] - updateterm); //<-- residue
	}
	
	if(buf.press_l1[i]<0)
		buf.press_l1[i]=0;

	float rhoil = buf.rho_adv[i]+updateterm+buf.aii[i]*buf.press_l[i];
	
	buf.densityResidue[i] = -buf.calcBuffer[i].restdens+(buf.rho_adv[i]+updateterm+buf.aii[i]*buf.press_l[i]);
	
	cfloat3 force(0,0,0);
	force = buf.dii[i]*buf.press_l[i] + buf.dijpj[i];
	//if (i%100==0) {
	//	printf("%f %f %f %f\n",force.x,force.y,force.z,buf.press_l[i]);
	//}
	buf.calcBuffer[i].accel = force/paramCarrier.dt /paramCarrier.dt;
	buf.calcBuffer[i].pressure = buf.press_l1[i];
}

__global__ void IntegrateIISPH(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;
	if (i >= pnum) return;
	if (buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;

	buf.calcBuffer[i].vel = buf.calcBuffer[i].veleval + buf.calcBuffer[i].accel * paramCarrier.dt;
	
	buf.displayBuffer[i].pos += buf.calcBuffer[i].vel * paramCarrier.dt/paramCarrier.simscale;
	//buf.calcBuffer[i].vel.Set(0,0,0);
}

__global__ void ComputeDensityIISPH(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if (i >= pnum) return;

	if (buf.displayBuffer[i].type==TYPE_BOUNDARY)
		return;

	// Get search cell
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range


	float sum = 0.0;
	float dens;

	dens = buf.calcBuffer[i].restdens;

	//Get Fluid Density
	for (int c=0; c < paramCarrier.neighbornum; c++) {
		contributeDensity(i, sum, gc + paramCarrier.neighborid[c], buf);
	}

	sum = sum * paramCarrier.kpoly6;

	if (sum == 0.0)
		sum = 1.0;

	buf.calcBuffer[i].dens = 1/sum;
}