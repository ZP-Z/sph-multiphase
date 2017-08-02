#include <assert.h>
#include <stdio.h>
#include <conio.h>


#include "fluid_system.h"
#include "fluid_system_host.cuh"



void FluidSystem::SetupAddBubble(cfloat3 min, cfloat3 max, float spacing, int constitution){
   /* cfloat3 pos;
	int n = 0, p;
	float dx, dy, dz, x, y, z;
	int cntx, cnty, cntz;
	cntx = ceil( (max.x-min.x) / spacing );
	cntz = ceil( (max.z-min.z) / spacing );
	int cnt = cntx * cntz;
	int xp, yp, zp, c2;
	float odd;

	if(constitution >= fluidnum)
		return;
	
	dx = max.x-min.x;
	dy = max.y-min.y;
	dz = max.z-min.z;
	float randx[3];
	float ranlen = 0.2;

	c2 = cnt/2;
	for (float y = min.y; y <= max.y; y += spacing ) {	
		for (int xz=0; xz < cnt; xz++ ) {
			randx[0] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
			randx[1] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
			randx[2] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;

			x = min.x + (xz % int(cntx))*spacing;
			z = min.z + (xz / int(cntx))*spacing;

			p = AddParticle ();
			if ( p != -1 ) {
				n++;
				(mPos+p)->Set ( x+randx[0],y+randx[1],z+randx[2]);

				*(m_alpha+p*MAX_FLUIDNUM+constitution) = 1.0f;
                *(m_alpha_pre+p*MAX_FLUIDNUM+constitution) = 1.0f;
				*(m_restMass+p) = hostCarrier.massArr[constitution];
				*(m_restDensity+p) = hostCarrier.densArr[constitution];
                mClr[p] = COLORA(m_alpha[p*MAX_FLUIDNUM],m_alpha[p*MAX_FLUIDNUM+1],m_alpha[p*MAX_FLUIDNUM+2],0.8);
				*(m_visc+p) = hostCarrier.viscArr[constitution];
				*(MF_type+p) = 2; 
			}
		}
	}	
	printf("%d fluid has %d particles\n",constitution,n);*/
}

//void FluidSystem::SetupAddBound(BI2Reader bi2reader,int boundtype)
//{
//	printf("%f %f %f %f\n",m_Param[PBMASS],m_Param[PBRESTDENSITY],m_Param[PBVISC],m_Param[PBSTIFF]);
//	float x,y,z;
//	for (int i = 0;i<bi2reader.info.nbound;i++)
//	{
//		x = bi2reader.info.Pos[i].x/m_Param[PSIMSCALE];
//		y = bi2reader.info.Pos[i].y/m_Param[PSIMSCALE];
//		z = bi2reader.info.Pos[i].z/m_Param[PSIMSCALE];
//		if( x < m_Vec[PVOLMIN].x || x>m_Vec[PVOLMAX].x || y<m_Vec[PVOLMIN].y || y> m_Vec[PVOLMAX].y || z<m_Vec[PVOLMIN].z || z>m_Vec[PVOLMAX].z)
//			continue;
//		//if (x<0)
//		//	continue;
//		int p = AddParticle();
//		if (p!=-1)
//		{
//			
//			(mPos+p)->Set (x,z,y);
//			*(mClr+p) = COLORA(1,1,1,1); 
//			*(mIsBound+p) = boundtype;
//			*(m_restMass+p) = m_Param[PBMASS];
//			*(m_restDensity+p) = m_Param[PBRESTDENSITY];
//			*(m_visc+p) = m_Param[PBVISC];
//			*(MF_type+p) = 2;//which means bubble (project-u)
//		}
//	}
//}
//
//void FluidSystem::SetupAddShape(BI2Reader bi2reader,int cat)
//{
//	printf("%f %f %f %f\n",m_Param[PBMASS],m_Param[PBRESTDENSITY],m_Param[PBVISC],m_Param[PBSTIFF]);
//	for (int i = 0;i<bi2reader.info.np;i++)
//	{
//		int p = AddParticle();
//		if (p!=-1)
//		{
//		//	printf("%f %f %f\n",bi2reader.info.Pos[i].x/m_Param[PSIMSCALE],bi2reader.info.Pos[i].y/m_Param[PSIMSCALE],bi2reader.info.Pos[i].z/m_Param[PSIMSCALE]);
//			(mPos+p)->Set (bi2reader.info.Pos[i].x/m_Param[PSIMSCALE],bi2reader.info.Pos[i].y/m_Param[PSIMSCALE],bi2reader.info.Pos[i].z/m_Param[PSIMSCALE]);
//			*(mClr+p) = COLORA(1,1,1,1); 
//			*(m_alpha+p*MAX_FLUIDNUM+cat) = 1.0f;*(m_alpha_pre+p*MAX_FLUIDNUM+cat) = 1.0f;
//			*(m_restMass+p) = hostCarrier.massArr[cat];
//			*(m_restDensity+p) = hostCarrier.densArr[cat];
//			*(m_visc+p) = hostCarrier.viscArr[cat];
//			*(mIsBound+p) = 2;
//			*(MF_type+p) = 1; //which means deformable (project-u)
//		}
//	}
//}

//void FluidSystem::SetupMfAddDeformVolume ( cfloat3 min, cfloat3 max, float spacing, cfloat3 offs, int type )
//{
//	cfloat3 pos;
//	int n = 0, p;
//	float dx, dy, dz, x, y, z;
//	int cntx, cnty, cntz;
//	cntx = ceil( (max.x-min.x-offs.x) / spacing );
//	cntz = ceil( (max.z-min.z-offs.z) / spacing );
//	cfloat3 mid = (min + max);
//	mid.x /= 2;
//	mid.y /= 2;
//	mid.z /= 2;
//	float radius = 0.5*(max.x - min.x);
//	int cnt = cntx * cntz;
//	int xp, yp, zp, c2;
//	float odd;
//
//	dx = max.x-min.x;
//	dy = max.y-min.y;
//	dz = max.z-min.z;
//	
//	c2 = cnt/2;
//
//	float xcenter = max.x - dx/2;
//	float ycenter = max.y - dy/2;
//	float omega = 0.0;
//	float rx,ry;
//	
//	for (float y = min.y+offs.y; y <= max.y; y += spacing ) {	
//		if( type==1 && y + spacing > max.y )
//			continue;
//		for (int xz=0; xz < cnt; xz++ ) {
//			x = min.x+offs.x + (xz % int(cntx))*spacing;
//			z = min.z+offs.z + (xz / int(cntx))*spacing;
//			
//			float dis = sqrtf((x - mid.x)*(x - mid.x) + (y - mid.y)*(y - mid.y) + (z - mid.z)*(z - mid.z));
//			if (dis > radius)
//				continue;
//			
//			p = AddParticle ();
//			if ( p != -1 ) {
//				n++;
//				(mPos+p)->Set ( x,y,z);
//                
//				//*(mClr+p) = COLORA( 0.25, +0.25 + (y-min.y)*.75/dy, 0.25 + (z-min.z)*.75/dz, 1);  // (x-min.x)/dx
//				*(m_alpha+p*MAX_FLUIDNUM+type) = 1.0f;
//                *(m_alpha_pre+p*MAX_FLUIDNUM+type) = 1.0f;
//                mClr[p] = COLORA(0.9,0.2,0.2,1);
//
//				*(m_restMass+p) = hostCarrier.massArr[0];
//				*(m_restDensity+p) = hostCarrier.densArr[0];			
//
//				*(m_visc+p) = hostCarrier.viscArr[0];
//				*(MF_type+p) = type; //1 means deformable
//
//				rx = x - xcenter;
//				ry = y - ycenter;
//				//mVel[p].Set( 0, -10, 0);
//				//mVelEval[p].Set( ry*omega, -rx*omega,0);
//			}
//		}
//	}	
//	printf("%d fluid has %d particles\n",0,n);
//}

void FluidSystem::SetupAddSolid( cfloat3 min, cfloat3 max, float spacing, int composition){
	
	//cfloat3 pos;
	//int n = 0, p;
	//float dx, dy, dz, x, y, z;
	//int cntx, cnty, cntz;
	//cntx = ceil( (max.x-min.x) / spacing );
	//cntz = ceil( (max.z-min.z) / spacing );
	//
	////rotation initialization
	//cfloat3 mid = (min + max);
	//mid.x /= 2;
	//mid.y /= 2;
	//mid.z /= 2;
	//float radius = 0.5*(max.x - min.x);
	//int cnt = cntx * cntz;
	//int xp, yp, zp, c2;
	//float odd;

	//dx = max.x-min.x;
	//dy = max.y-min.y;
	//dz = max.z-min.z;
	//
	//c2 = cnt/2;

	//float xcenter = max.x - dx/2;
	//float ycenter = max.y - dy/2;
	//float zcenter = max.z - dz/2;
	//float omega = 1.0;
	//float rx,ry;

	//float randx[3]={0,0,0};
	//float ranlen = 0.2;
	//int xi, yi, zi;

	//for (float y = min.y; y <= max.y; y += spacing ) {	

	//	for (int xz=0; xz < cnt; xz++ ) {
	//		xi = xz % cntx;
	//		zi = xz / cntx;
	//		yi = y / spacing;

	//		x = min.x + (xz % int(cntx))*spacing;
	//		z = min.z + (xz / int(cntx))*spacing;

	//		//if (((x-min.x)*(x-min.x)+(y-ycenter)*(y-ycenter)+(z-min.z)*(z-min.z))<225)
	//		//	continue;

	//		
	//		
	//		/*float dis = sqrtf((x - mid.x)*(x - mid.x) + (y - mid.y)*(y - mid.y) + (z - mid.z)*(z - mid.z));
	//		if (dis > radius)
	//			continue;*/
	//		
	//		p = AddParticle ();
	//		if ( p != -1 ) {
	//			n++;

	//			//randx[0] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
	//			//randx[1] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
	//			//randx[2] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;

	//			(mPos+p)->Set ( x+randx[0],y+randx[1],z+randx[2]);
 //               
	//			*(m_alpha     + p*MAX_FLUIDNUM + composition) = 1.0f;
 //               *(m_alpha_pre + p*MAX_FLUIDNUM + composition) = 1.0f;
 //               mClr[p] = COLORA( m_alpha[p*MAX_FLUIDNUM], m_alpha[p*MAX_FLUIDNUM+1], m_alpha[p*MAX_FLUIDNUM+2],1);

	//			m_restMass[p] = hostCarrier.massArr[composition];
	//			m_restDensity[p] = hostCarrier.densArr[composition];
	//			m_visc[p] = hostCarrier.viscArr[composition];
	//			MF_type[p] = 1; //1 means deformable

	//			//rx = x - xcenter;
	//			//ry = y - ycenter;
	//			//mVel[p].Set( 0, -1, 0);
	//			//mVel[p].Set( 0, , 0);
	//		}
	//	}
	//}	
	//printf("%d solid has %d particles\n",0,n);
}

void read(FILE* fp, cfloat3& dst){
    fscanf(fp, "%f,%f,%f", &dst.x, &dst.y, &dst.z);
    return;
}
void read(FILE* fp, float& dst) {
	fscanf(fp, "%f", &dst);
	return;
}

void FluidSystem::CalcBoundarySpacing() {
	float bmass = hostCarrier.bmass;
	float bdens = hostCarrier.bRestdensity;

	boundarySpacing = pow((bmass / bdens),1/3.0f)/hostCarrier.simscale;
	printf("boundaryspacing %f\n",boundarySpacing);
}

void FluidSystem::LoadBoundary(std::string boundfile){

	CalcBoundarySpacing();

    FILE* fp = fopen(boundfile.c_str(), "r");
    if (fp==NULL){
        printf("null boundary file");
		return;
	}
    boundaryPArray.clear();
    
    char buffer[1000];
    cfloat3 min,max;
    while (true){
        if(fscanf(fp, "%s", buffer)==EOF)
            break;
        if (!strcmp(buffer, "WALL")){
            read(fp, min);
            read(fp, max);
            SetupAddWall(min,max);
        }
		if (!strcmp(buffer, "OPEN_BOX")) {
			read(fp,min);
			read(fp,max);
			float thickness;
			read(fp, thickness);
			SetupAddOpenBox(min,max,thickness);
		}
    }
    return;
}

void FluidSystem::SetupAddWall(cfloat3 min, cfloat3 max){
    cfloat3 pos;
	int n = 0, p;
	float dx, dy, dz, x, y, z;
	int cntx, cnty, cntz;

    float spacing = boundarySpacing;
	cntx = ceil( (max.x-min.x) / spacing );
	cntz = ceil( (max.z-min.z) / spacing );
	
	//rotation initialization
	
	int cnt = cntx * cntz;
	int xp, yp, zp, c2;
	float odd;

	dx = max.x-min.x;
	dy = max.y-min.y;
	dz = max.z-min.z;
	
	c2 = cnt/2;

	float randx[3]={0,0,0};
	float ranlen = 0.2;

	for (float y = min.y; y <= max.y; y += spacing ) {	
		for (int xz=0; xz < cnt; xz++ ) {
			x = min.x + (xz % int(cntx))*spacing;
			z = min.z + (xz / int(cntx))*spacing;
			
			p = AddParticle ();
			if ( p != -1 ) {
				n++;
				/*
				randx[0] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
				randx[1] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
				randx[2] = (float)rand()/RAND_MAX * ranlen - ranlen*0.5;
				*/
				displayBuffer[p].pos.Set( x,y,z);
				if(z>0)
					displayBuffer[p].color.Set(1,1,1,0);
				else
					displayBuffer[p].color.Set(1, 1, 1, 0.1);
				displayBuffer[p].type = 1;

				calculationBuffer[p].mass = hostCarrier.bmass;
				calculationBuffer[p].restdens = hostCarrier.bRestdensity;
				calculationBuffer[p].visc = hostCarrier.bvisc;
			}
		}
	}	
}

void FluidSystem::SetupAddOpenBox(cfloat3 min, cfloat3 max, float thickness) {
	//ground
	SetupAddWall(cfloat3(min.x,min.y-thickness,min.z), cfloat3(max.x, min.y, max.z));
	//front
	SetupAddWall(cfloat3(min.x-thickness,min.y,min.z-thickness), cfloat3(max.x+thickness, max.y, min.z));
	//back
	SetupAddWall(cfloat3(min.x-thickness, min.y, max.z), cfloat3(max.x+thickness, max.y, max.z+thickness));
	//left
	SetupAddWall(cfloat3(min.x-thickness, min.y, min.z), cfloat3(min.x, max.y, max.z));
	//right
	SetupAddWall(cfloat3(max.x, min.y, min.z), cfloat3(max.x+thickness, max.y, max.z));
}

void FluidSystem::SetupAddSphere(cfloat3 center, float radius){

}

void FluidSystem::SetupAddShell(cfloat3 center, float innner, float outter, float cutlevel){

}
