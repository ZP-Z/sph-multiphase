

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

__device__ FluidParams	simData;
__device__ ParamCarrier paramCarrier;


void CarryParam(ParamCarrier& hostCarrier){
	cudaMemcpyToSymbol( paramCarrier, &hostCarrier, sizeof(hostCarrier));
}

void updateParam(FluidParams* paramCPU){
    cudaMemcpyToSymbol ( simData, paramCPU, sizeof(FluidParams) );
}

__device__ void findNearest ( int i, float& mindis, int cell, bufList buf )
{			
	cfloat3 dist;
	float dsq;
	cfloat3 p = buf.mpos[i];

	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;
	int j;

	if ( buf.mgridcnt[cell] == 0 ) return ;
	
	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];
	for ( int cndx = cfirst; cndx < clast; cndx++ ) {
		j = buf.mgrid[cndx];
		
		if (buf.misbound[j] == 0)
		{
			j = buf.mgrid[cndx];
			dist = p - buf.mpos[ buf.mgrid[cndx] ];
			dsq = dot(dist,dist);

			if ( dsq < r2 && dsq > 0.0 && dsq*d2<mindis) 
			{
				mindis = dsq*d2;
				buf.midsort[i] = j;
			} 
		}
	}
	
	return ;
}

//Sorting

//get particle index id
__global__ void InitialSort ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;
	register cfloat3 gridMin = simData.gridMin;
	register cfloat3 gridDelta = simData.gridDelta;
	register cint3 gridRes = simData.gridRes;
	register cint3 gridScan = simData.gridScanMax;
//	register float poff = simData.psmoothradius / simData.psimscale;

	register int		gs;
	register cfloat3		gcf;
	register cint3		gc;

	gcf = (buf.mpos[i] - gridMin) * gridDelta; 
	gc = cint3( int(gcf.x), int(gcf.y), int(gcf.z) );
	/*if(buf.MFid[i]==10000){
		printf("%f %f %f %f %f %f\n",buf.mpos[i].x,buf.mpos[i].y,buf.mpos[i].z,
			gridMin.x,gridMin.y,gridMin.z);
		printf("gc %d %d %d\n",gc.x,gc.y,gc.z);
	}*/
	gs = (gc.y * gridRes.z + gc.z)*gridRes.x + gc.x;
	if ( gc.x >= 1 && gc.x <= gridScan.x && gc.y >= 1 && gc.y <= gridScan.y && gc.z >= 1 && gc.z <= gridScan.z ) {
		buf.mgcell[i] = gs;											// Grid cell insert.
		buf.midsort[i] = i;
//		buf.mgndx[i] = atomicAdd ( &buf.mgridcnt[ gs ], 1 );		// Grid counts.
//		gcf = (-cfloat3(poff,poff,poff) + buf.mpos[i] - gridMin) * gridDelta;
//		gc = make_int3( int(gcf.x), int(gcf.y), int(gcf.z) );
//		gs = ( gc.y * gridRes.z + gc.z)*gridRes.x + gc.x;
		//buf.mcluster[i] = gs;				-- make sure it is allocated!
	} else {
		buf.mgcell[i] = GRID_UNDEF;
		buf.midsort[i] = i;
		//buf.mcluster[i] = GRID_UNDEF;		-- make sure it is allocated!
	}
	/*if (buf.MFid[i]==0) {
		printf("i %d\n", buf.mgcell[i]);
	}*/
}

//markup the head and tail of each cell
__global__ void CalcFirstCnt ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;
	if (i>=pnum) return;
	if ((i == 0 || buf.mgcell[i]!=buf.mgcell[i-1]))
	{
		if (buf.mgcell[i]!=GRID_UNDEF)buf.mgridoff[buf.mgcell[i]] = i;
	}
	__syncthreads();
	if (i!=0 && buf.mgcell[i]!=buf.mgcell[i-1] && buf.mgcell[i-1]!=GRID_UNDEF)
		buf.mgridcnt[buf.mgcell[i-1]] = i;
	if (i == pnum-1 && buf.mgcell[i]!=GRID_UNDEF)
		buf.mgridcnt[buf.mgcell[i]] = i + 1;
	
}
__global__ void GetCnt ( bufList buf, int pnum )
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;		// particle index
	if (i>=pnum) return ;
	if (buf.mgcell[i]!=GRID_UNDEF)
	{
		buf.mgndx[i] = i - buf.mgridoff[buf.mgcell[i]];
		if (buf.mgndx[i] == 0)
			buf.mgridcnt[buf.mgcell[i]] -= buf.mgridoff[buf.mgcell[i]];
	}
}

//deep copy sort
__global__ void CountingSortFull_ ( bufList buf, int pnum )
{
	//for each new position, find old particle and read value
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;		// particle index				
	if ( i >= pnum ) return;

	//for each old particle, find new position and write
	//uint icell = *(uint*) (buf.msortbuf + pnum*BUF_GCELL + i*sizeof(uint) );
	//uint indx =  *(uint*) (buf.msortbuf + pnum*BUF_GNDX + i*sizeof(uint) );
	//int sort_ndx = buf.mgridoff[ icell ] + indx;				// global_ndx = grid_cell_offet + particle_offset
	//i = buf.midsort[i];
	
	int exId = buf.midsort[i];
	int cell = buf.mgcell[i];
	int sort_ndx = i;
	i = exId;

	
	if ( cell != GRID_UNDEF ) {
		buf.mgrid[ sort_ndx ] = sort_ndx;			// full sort, grid indexing becomes identity
		char* bpos = buf.msortbuf + i*sizeof(cfloat3);
		buf.mpos[ sort_ndx ] =		*(cfloat3*) (bpos);
		buf.mvel[ sort_ndx ] =		*(cfloat3*) (bpos + pnum*BUF_VEL );
		buf.mveleval[ sort_ndx ] =	*(cfloat3*) (bpos + pnum*BUF_VELEVAL );
		buf.mforce[sort_ndx] = *(cfloat3*)(bpos + pnum*BUF_FORCE);
		buf.mpress[sort_ndx] = *(float*)(buf.msortbuf + pnum*BUF_PRESS + i*sizeof(float));
		//buf.last_mpress[sort_ndx] = *(float*)(buf.msortbuf + pnum*BUF_LAST_PRESS + i*sizeof(float));
		buf.mdensity[ sort_ndx ] =	*(float*) (buf.msortbuf + pnum*BUF_DENS + i*sizeof(float) );
		buf.mclr[ sort_ndx ] =		*(uint*) (buf.msortbuf + pnum*BUF_CLR+ i*sizeof(uint) );		// ((uint) 255)<<24; -- dark matter
		buf.misbound[ sort_ndx ] =		*(int*) (buf.msortbuf + pnum*BUF_ISBOUND+ i*sizeof(int) );		// ((uint) 255)<<24; -- dark matter
		buf.accel[sort_ndx] = *(cfloat3*)(buf.msortbuf + pnum*BUF_ACCEL + i*sizeof(cfloat3));

		buf.mgcell[ sort_ndx ] =	cell;

		//multi fluid
		int mul_sort_ndx = sort_ndx*MAX_FLUIDNUM;
		for( uint fcount = 0; fcount < simData.mf_catnum; fcount++)
		{
			//char* bmul = buf.msortbuf + i*sizeof(float)*MAX_FLUIDNUM + fcount * sizeof(float);
			buf.mf_alpha[mul_sort_ndx+fcount] =			*(float*)(buf.msortbuf +  pnum*BUF_ALPHA +   i*sizeof(float)*MAX_FLUIDNUM + fcount * sizeof(float));
			buf.mf_alpha_pre[mul_sort_ndx+fcount] =		*(float*)(buf.msortbuf +  pnum*BUF_ALPHAPRE+ i*sizeof(float)*MAX_FLUIDNUM + fcount * sizeof(float));
			//buf.mf_pressure_modify[mul_sort_ndx+fcount] =	*(float*)(bmul + pnum*BUF_PRESSMODI);
			buf.mf_vel_phrel[mul_sort_ndx+fcount] =		*(cfloat3*)(buf.msortbuf + pnum*BUF_VELPHREL +  i*sizeof(cfloat3)*MAX_FLUIDNUM + fcount*sizeof(cfloat3));
			buf.mf_alphagrad[mul_sort_ndx+fcount] =		*(cfloat3*)(buf.msortbuf + pnum*BUF_ALPHAGRAD + i*sizeof(cfloat3)*MAX_FLUIDNUM + fcount*sizeof(cfloat3));
		}
		buf.mf_pressure_modify[ sort_ndx ] = *(float*) (buf.msortbuf + pnum*BUF_PRESSMODI + i*sizeof(float));
		buf.mf_restmass[ sort_ndx ] = *(float*) (buf.msortbuf + pnum*BUF_RMASS + i*sizeof(float));
		buf.mf_restdensity[sort_ndx] = *(float*) (buf.msortbuf + pnum*BUF_RDENS + i*sizeof(float));
		buf.mf_visc[sort_ndx] = *(float*) (buf.msortbuf + pnum*BUF_VISC + i*sizeof(float));
		buf.mf_velxcor[sort_ndx] = *(cfloat3*)(buf.msortbuf + pnum*BUF_VELXCOR + i*sizeof(cfloat3));
		buf.MFtype[sort_ndx] = *(int*)(buf.msortbuf+ pnum*BUF_INDICATOR + i*sizeof(int));
		for(int k=0; k<9; k++){
			buf.MFtensor[sort_ndx*9+ k] = *(float*)(buf.msortbuf + pnum*BUF_TENSOR + i*sizeof(float)*9 + k*sizeof(float));
			buf.MFtemptensor[sort_ndx*9+k] = *(float*)(buf.msortbuf + pnum*BUF_TEMPTENSOR + i*sizeof(float)*9 + k*sizeof(float));
			buf.MFRtensor[sort_ndx*9+k] = *(float*)(buf.msortbuf + pnum*BUF_RTENSOR + i*sizeof(float)*9 + k*sizeof(float));
		}
		buf.MFid[sort_ndx] = *(int*)(buf.msortbuf + pnum*BUF_BORNID + i*sizeof(int));
		buf.MFidTable[ buf.MFid[sort_ndx]] = sort_ndx;
	}
	else{
		buf.mgcell[sort_ndx] = GRID_UNDEF;
		buf.mpos[sort_ndx] = cfloat3(-1000, -1000, -1000);
	}
}


__global__ void mfFindNearest (bufList buf,int pnum)
{
	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;				
	if ( i >= pnum ) return;
	
	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;
	gc -= nadj;
	
	// Sum Pressures
	cfloat3 pos = buf.mpos[ i ];
	float mindis = 65535;

	if (buf.misbound[i]==1) //boundary particles
	{
		buf.midsort[i] = i;
		buf.mf_restmass[i] = simData.pmass;
		for (int c = 0; c<simData.gridAdjCnt; c++)
		{
			findNearest(i,mindis,gc+simData.gridAdj[c],buf);
		}
		if (buf.midsort[i]!=i)
			buf.mf_restmass[i] = buf.mf_restmass[buf.midsort[i]];

	}
}

__device__ void contributeDensity_boundary(uint i, float& res, uint cell, bufList buf) {
	cfloat3 dist;
	cfloat3 p = buf.mpos[i];
	float dsq, c;
	float massj;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	int j;

	if (buf.mgridcnt[cell] == 0)
		return;

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];
		if ( buf.misbound[j] != 1 )
			continue;

		dist = p - buf.mpos[j];

		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if (dsq < r2 && dsq > 0.0) {
			c = (r2 - dsq)*d2;
			res += c * c * c ;//* buf.mf_restmass[j];
		}
	}
	return;
}

__device__ void contributeDensity (uint i, float& res, uint cell, bufList buf){
	cfloat3 dist;
	cfloat3 p = buf.mpos[i];
	float dsq, c;
	float massj;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	int j;

	if ( buf.mgridcnt[cell] == 0 )
		return;
	/*if(buf.MFid[i]==0){
		printf("%d %d\n", cell, buf.mgridcnt[cell]); 
	}*/
	int cfirst = buf.mgridoff[ cell ];
	int clast = cfirst + buf.mgridcnt[ cell ];

	for ( int cndx = cfirst; cndx < clast; cndx++ ) {
		j = buf.mgrid[cndx];
        //if ( buf.misbound[j]==0)
		massj = buf.mf_restmass[j];
		//else
		//	massj = buf.mf_restdensity[i] * buf.mdensity[j] * simData.fsa; // fluid density * boundary volume

		dist = p - buf.mpos[j];

		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		if ( dsq < r2 && dsq > 0.0) {
			c = (r2 - dsq)*d2;
			res += c * c * c * massj;	
		} 
	}
	return;
}

__global__ void ComputeBoundaryVolume(bufList buf, int pnum) {
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if (i >= pnum) return;

	if(buf.misbound[i]==0)
		return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[i];
	if (gc == GRID_UNDEF) return;						// particle out-of-range
	gc -= nadj;

	float sum = 0.0;

	//Get Boundary Density
	
	for (int c=0; c < simData.gridAdjCnt; c++) {
		contributeDensity_boundary(i, sum, gc + simData.gridAdj[c], buf);
	}

	sum += simData.r2 * simData.r2 * simData.r2 ;//* buf.mf_restmass[i];
	sum = sum * simData.poly6kern;

	if (sum == 0.0) {
		printf("boundary density zero error.\n");
		sum = 1.0;
	}
		
	buf.mdensity[i] = 1 / sum ; //actually the volume
}

__global__ void ComputeDensityPressure(bufList buf,int pnum){
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum ) return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	float sum = 0.0;
	float dens;
	if(buf.misbound[i]==1)
		dens = buf.mf_restdensity[ buf.midsort[i] ];
	else
		dens = buf.mf_restdensity[i];

	//Get Fluid Density
	if(buf.misbound[i]==0)
		for (int c=0; c < simData.gridAdjCnt; c++) {
			contributeDensity(i, sum, gc + simData.gridAdj[c], buf);
		}
	
	/*if (buf.midsort[i]==0) {
		printf("contribute? %f %f\n", sum, dens);
	}*/

	sum += simData.r2 * simData.r2 * simData.r2 * buf.mf_restmass[i];
	sum = sum * simData.poly6kern;
	
	if ( sum == 0.0 )
		sum = 1.0;
	
	buf.mdensity[i] = 1/sum;

	//buf.mpress[ i ] = ( sum - dens ) * simData.pintstiff;
	//buf.mpress[ i ] = (pow( sum/dens,7.0f )-1) * simData.pintstiff;
    buf.mpress[i] = 2.5 * dens * (pow(sum / dens, 7.0f) - 1);

    if (buf.mpress[i]<0)
        buf.mpress[i] = 0;
	//if (buf.mpress[i]>65535)
	//	printf("press i too big error");

	buf.mclr[i] = COLORA(1, 1-buf.mpress[i]/1000, 1-buf.mpress[i]/1000, 1);
}






__device__ cfloat3 contributeForce_new(int i, int muli, cfloat3 ipos, cfloat3 iveleval, float ipress, float idens, int cell, bufList buf, cfloat3* ivelxcor, float ivisc)
{
	//Force here represents the acceleration
	float dsq, c;
	register float d2 = simData.psimscale * simData.psimscale;
	register float r2 = simData.r2/d2;

	cfloat3 dist, vmr;
	float cmterm, cmterm1;

	float pmterm, vmterm;

	int j;
	float aveDenij, cx, xterm;

	if (buf.mgridcnt[cell] == 0) return cfloat3(0, 0, 0);

	int cfirst = buf.mgridoff[cell];
	int clast = cfirst + buf.mgridcnt[cell];

	cfloat3 force = cfloat3(0, 0, 0);

	for (int cndx = cfirst; cndx < clast; cndx++) {
		j = buf.mgrid[cndx];
		//if (buf.misbound[j]==1)
		//    continue;

		dist = (ipos - buf.mpos[j]);		
		dsq = (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z);
		dist *= simData.psimscale;

		if (!(dsq < r2 && dsq > 0))
			continue;

		cx = (r2-dsq)*d2;
		dsq = sqrt(dsq*d2);
		c = (simData.psmoothradius - dsq);

		cmterm1 = simData.spikykern * c * c / dsq;
		cmterm = cmterm1 * buf.mf_restmass[j] * buf.mdensity[j];


		if (buf.misbound[j] != 1) //fluid & solid, fluid & fluid, solid & solid
		{
			//pressure
			pmterm = - 0.5* cmterm * (buf.mpress[j]) *idens;
			force += dist * pmterm;

			//viscosity
			vmr = iveleval - buf.mveleval[j];
			vmterm = cmterm * (ivisc+buf.mf_visc[j]) * idens;
			force += vmr*vmterm;

		}
			
		//else { //boundary & non-boundary
		//	//pressure
		//	pmterm = - cmterm1 * buf.mf_restdensity[i] * buf.mdensity[j] *  (ipress ) *idens*idens;
		//	force += pmterm * dist;
		//		
		//	vmr = iveleval - buf.mveleval[j];
		//	vmterm = cmterm * (ivisc+buf.mf_visc[j]) * idens;
		//	force += vmterm * vmr;
		//		
		//	//artificial boundary viscosity			
		//	vmr = iveleval - buf.mveleval[j];
		//	float pi_ij = vmr.x*dist.x + vmr.y*dist.y + vmr.z*dist.z;
		//	if (pi_ij < 0) {
		//		pi_ij = pi_ij / (dist.x*dist.x + dist.y*dist.y + dist.z*dist.z + r2 * 0.01);
		//		pi_ij = pi_ij * 2 * simData.psmoothradius * (ivisc + buf.mf_visc[j]) * idens /2;
		//		pi_ij = - cmterm1 * buf.mf_restdensity[i] * buf.mdensity[j] * pi_ij;
		//		force += (pi_ij * dist * simData.visc_factor);
		//	}
		//}
		
	}
	return force;
}

__global__ void ComputeForce (bufList buf, int pnum){
	uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum)
		return;
	if(buf.misbound[i]==1)
		return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	register uint muloffseti = i * MAX_FLUIDNUM;
	register cfloat3 ipos = buf.mpos[ i ];
	register cfloat3 iveleval = buf.mveleval[ i ];
	register float ipress = buf.mpress[ i ];
	register float idens = buf.mdensity[ i ];
	register cfloat3 *ivelxcor = buf.mf_velxcor+i;
	register float ivisc = buf.mf_visc[i];

	register cfloat3 force = cfloat3(0,0,0);	
	*ivelxcor = cfloat3(0,0,0);

	for (int c=0; c < simData.gridAdjCnt; c++) {
		force += contributeForce_new ( i, muloffseti, ipos, iveleval, ipress, idens, gc + simData.gridAdj[c], buf, ivelxcor, ivisc);
	}
	buf.mforce[ i ] = force;
	if (buf.MFid[i]==0) {

		//buf.mclr[i] = COLORA(0,0,1,1);
		//printf("%f\n",buf.mpress[i]);
		//printf("%f %f %f\n",buf.mforce[i].x, buf.mforce[i].y, buf.mforce[i].z);

	}
}

__device__ void contributeCover(int i, int cell, bufList buf){

}

__global__ void SurfaceDetection(bufList buf, int pnum){
    uint i = blockIdx.x * blockDim.x + threadIdx.x;	// particle index				
	if ( i >= pnum)
		return;
	if(buf.misbound[i]==1)
		return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	register uint muloffseti = i * MAX_FLUIDNUM;
	register cfloat3 ipos = buf.mpos[ i ];
	register cfloat3 iveleval = buf.mveleval[ i ];
	register float ipress = buf.mpress[ i ];
	register float idens = buf.mdensity[ i ];
	register cfloat3 *ivelxcor = buf.mf_velxcor+i;
	register float ivisc = buf.mf_visc[i];

	register cfloat3 force = cfloat3(0,0,0);	
	*ivelxcor = cfloat3(0,0,0);

	for (int c=0; c < simData.gridAdjCnt; c++) {
		contributeCover ( i, gc + simData.gridAdj[c], buf );
	}
}


//Yang Tao Method
__device__ void contributeSurfaceTensionYT(int i, cfloat3& res, int cell, bufList buf)
{
    //Force here represents the acceleration
	float dsq, c;
	float simscale2 = simData.psimscale * simData.psimscale;	
	float sr2 = simData.r2 / simscale2;
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

	for ( int cndx = cfirst; cndx < clast; cndx++ ) {	
		j = buf.mgrid[ cndx ];	
        if (buf.misbound[j]==1)
            continue;

        dist = buf.mpos[i] - buf.mpos[j];
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
	if(buf.misbound[i]==1)
		return;

	// Get search cell
	int nadj = (1*simData.gridRes.z + 1)*simData.gridRes.x + 1;
	uint gc = buf.mgcell[ i ];
	if ( gc == GRID_UNDEF ) return;						// particle out-of-range
	gc -= nadj;

	register uint muloffseti = i * MAX_FLUIDNUM;
	register cfloat3 ipos = buf.mpos[ i ];
	register cfloat3 iveleval = buf.mveleval[ i ];
	register float ipress = buf.mpress[ i ];
	register float idens = buf.mdensity[ i ];
	register cfloat3 *ivelxcor = buf.mf_velxcor+i;
	register float ivisc = buf.mf_visc[i];

	register cfloat3 force = cfloat3(0,0,0);	

	for (int c=0; c < simData.gridAdjCnt; c++) {
		contributeSurfaceTensionYT ( i, force, gc + simData.gridAdj[c], buf);
	}

	buf.mforce[ i ] = force + buf.mforce[i] * paramCarrier.fBlending;

}

__global__ void AdvanceParticles(float time, float dt, float ss, bufList buf, int numPnts)
{

	uint i = __mul24(blockIdx.x, blockDim.x) + threadIdx.x;	// particle index				
	if (i >= numPnts) return;

	if (buf.mgcell[i] == GRID_UNDEF) {
		buf.mpos[i] = cfloat3(-1000, -1000, -1000);
		buf.mvel[i] = cfloat3(0, 0, 0);
		return;
	}

	// Get particle vars
	register cfloat3 accel, norm;
	register float diff, adj, accelLen;
	register cfloat3 pos = buf.mpos[i];
	register cfloat3 veval = buf.mveleval[i];
	
	accel = cfloat3(0, 0, 0);

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

		/*
		float frifac, press;
		press = buf.MFtensor[i*9]+buf.MFtensor[i*9+4]+buf.MFtensor[i*9+8];
		if (press < 0.0f)
			press = 0.0f;
		float xzmod = sqrt(veval.x*veval.x + veval.z*veval.z);
		if (xzmod > simData.sleepvel) {
			frifac = simData.bdamp * press;
			if (frifac < 15)
				frifac = 15;
			cfloat3 dv = cfloat3(-veval.x/xzmod, 0, -veval.z/xzmod) * frifac;
			accel += dv;
		}
		else {
			veval.x = 0;
			veval.z = 0;
		}*/
	
	}

	diff = simData.pradius - (simData.pboundmax.y - pos.y)*ss;

	if (diff > EPSILON) {
		norm = cfloat3(0, -1, 0);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;
	}

	// X-axis
	diff = simData.pradius - (pos.x - (simData.pboundmin.x + (sin(time*simData.pforce_freq)+1)*0.5 * simData.pforce_min))*ss;
	//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if (diff > EPSILON) {
		norm = cfloat3(1, 0, 0);
		adj = (simData.pforce_min+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}
	diff = simData.pradius - ((simData.pboundmax.x - (sin(time*simData.pforce_freq)+1)*0.5*simData.pforce_max) - pos.x)*ss;
	//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if (diff > EPSILON) {
		norm = cfloat3(-1, 0, 0);
		adj = (simData.pforce_max+1) * simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}

	// Z-axis
	diff = simData.pradius - (pos.z - simData.pboundmin.z) * ss;
	//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if (diff > EPSILON) {
		norm = cfloat3(0, 0, 1);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}
	diff = simData.pradius - (simData.pboundmax.z - pos.z)*ss;
	//	if (diff>simData.pradius) diff += simData.pradius*1000;
	if (diff > EPSILON) {
		norm = cfloat3(0, 0, -1);
		adj = simData.pextstiff * diff - simData.pdamp * dot(norm, veval);
		norm *= adj; accel += norm;//*scale_dens;
	}
	
	accel += buf.mforce[i];
	if(buf.MFid[i]==0){
		//printf("%f %f %f\n",buf.mforce[i].x, buf.mforce[i].y, buf.mforce[i].z);
		//printf("dens %f\n",buf.mpress[i]);
		//printf("pos %f %f %f\n",buf.mpos[i].x,buf.mpos[i].y,buf.mpos[i].z);
	}
	accel += simData.pgravity;
	// Accel Limit

	accelLen = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
	if ( accelLen > simData.AL2 ) {
	accel *= simData.AL / sqrt(accelLen);
	}
	cfloat3 vel = buf.mvel[i];
	if (buf.misbound[i] != 1)
	{

		//Velocity Limit
		/*float vmod = sqrtf(dot(vel, vel));
		if (vmod > simData.fluidVConstraint)
			vel *= simData.fluidVConstraint / vmod;*/

		cfloat3 vnext = accel*dt + vel;		// v(t+1/2) = v(t-1/2) + a(t) dt
		buf.mpos[i] += vnext * dt / ss;
		buf.mveleval[i] = (vel+vnext)*0.5;
		buf.mvel[i] = vnext;
		buf.accel[i] = accel;

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
}
