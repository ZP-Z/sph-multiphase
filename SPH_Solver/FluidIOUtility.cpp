
#include <assert.h>
#include <stdio.h>
#include <conio.h>

#include "fluid_system.h"
#include "fluid_system_host.cuh"

//void FluidSystem::outputepsilon(FILE* fp){
//	//fprintf(fp,"Frame: %d\n",nOutFrame);
//	for (int i = 0;i<NumPoints();i++)
//		fprintf(fp,"%d %f %f\n",i,MF_pepsilon[i],mPressure[i] + m_pressure_modify[i]);
//	fflush(stdout);
//	nOutFrame++;
//}

void FluidSystem::outputFile()
{
    char dsttmp[100];
	FILE* fp;
	sprintf(dsttmp,"OutputData\\data_%04d.txt",frameNo);
	fp = fopen(dsttmp,"w");
	fprintf(fp,"%d\n", pointNum);
	cfloat3* ppos = mPos;
	for (int i = 0;i<pointNum;i++,ppos++){
		if (mIsBound[i]!=1){
			fprintf(fp,"%f %f %f ",mPos[i].x,mPos[i].y,mPos[i].z); //x,y,z
			for(int k=0; k<MAX_FLUIDNUM; k++)
				fprintf(fp, "%f ", m_alpha[i*MAX_FLUIDNUM+k]); //alpha1, alpha2, alpha3
			fprintf(fp,"%f %f %d\n", m_restMass[i], m_restDensity[i],MF_type[i]); //mass, density, type
		}
	}
	fclose(fp);
	outputNo++;
}


//void FluidSystem::DrawGrid ()
//{
//	cfloat3 gd (1, 1, 1);
//	cfloat3 gc;
//	gd /= m_GridDelta;		
//	
//	glBegin ( GL_LINES );	
//	for (int z=0; z <= m_GridRes.z; z++ ) {
//		for (int y=0; y <= m_GridRes.y; y++ ) {
//			gc.Set ( 1, y, z);	gc /= m_GridDelta;	gc += m_GridMin;
//			glVertex3f ( m_GridMin.x, gc.y, gc.z );	glVertex3f ( m_GridMax.x, gc.y, gc.z );
//		}
//	}
//	for (int z=0; z <= m_GridRes.z; z++ ) {
//		for (int x=0; x <= m_GridRes.x; x++ ) {
//			gc.Set ( x, 1, z);	gc /= m_GridDelta;	gc += m_GridMin;
//			glVertex3f ( gc.x, m_GridMin.y, gc.z );	glVertex3f ( gc.x, m_GridMax.y, gc.z );
//		}
//	}
//	for (int y=0; y <= m_GridRes.y; y++ ) {
//		for (int x=0; x <= m_GridRes.x; x++ ) {
//			gc.Set ( x, y, 1);	gc /= m_GridDelta;	gc += m_GridMin;
//			glVertex3f ( gc.x, gc.y, m_GridMin.z );	glVertex3f ( gc.x, gc.y, m_GridMax.z );
//		}
//	}
//	glEnd ();
//}
//
//void FluidSystem::DrawText ()
//{
//	char msg[100];
//
//	
//	cfloat3* ppos = mPos;
//	DWORD* pclr = mClr;
//	cfloat3 clr;
//	for (int n = 0; n < NumPoints(); n++) {
//	
//		sprintf ( msg, "%d", n );
//		glColor4f ( (RED(*pclr)+1.0)*0.5, (GRN(*pclr)+1.0)*0.5, (BLUE(*pclr)+1.0)*0.5, ALPH(*pclr) );
//		drawText3D ( ppos->x, ppos->y, ppos->z, msg );
//		ppos++;
//		pclr++;
//	}
//}
//
//void FluidSystem::Draw ( Camera3D& cam, float rad )
//{
//	char* dat;
//	cfloat3* ppos;
//	float* pdens;
//	DWORD* pclr;
//		
//
//	glDisable ( GL_LIGHTING );
//
//	switch ( (int) m_Param[PDRAWGRID] ) {
//	case 0:
//		break;
//	case 1: 
//		glColor4f ( 0.7, 0.7, 0.7, 0.05 );
//		DrawGrid ();
//		break;
//	};
//
//	if ( m_Param[PDRAWTEXT] == 1.0 ) {
//		DrawText ();
//	};
//
//	// Draw Modes
//	// DRAW_POINTS		0
//	// DRAW_SPRITES		1
//	
//	switch ( (int) m_Param[PDRAWMODE] ) {
//	
//	case 1: //actually used
//
//
//		//Control Visability ==========================
//		//hide half bound
//		//for (int i = 0; i<NumPoints(); i++)
//		//	if (mIsBound[i]==1 /*&& mPos[i].x>0*/ )
//		//		mPos[i].x=mPos[i].y=mPos[i].z=-1000;
//		
//		switch(GetSimuP(SHOW_TYPE)){ //controlled by 7
//		case 1://bubble only
//			for (int i = 0;i<NumPoints();i++)
//				if (MF_type[i]!=2)
//					mPos[i].x=mPos[i].y=mPos[i].z=-1000;
//			break;
//		}
//
//
//		//Control Color   ==============================
//
//		switch (GetSimuP(SHOW_PROPERTY)) {
//		case 0://alpha
//			float r, g, b;
//			for (int i=0; i<NumPoints(); i++){
//				if (mIsBound[i]==0) {
//					r = 1 - (m_alpha[i*3+1]+m_alpha[i*3+2])*0.5;
//					g = 1 - (m_alpha[i*3]  +m_alpha[i*3+2])*0.5;
//					b = 1 - (m_alpha[i*3]  +m_alpha[i*3+1])*0.5;
//					mClr[i] = COLORA(r, g, b, 1);
//				}	
//				else
//					mClr[i] = COLORA(1,1,1,0.5);
//
//				//if(MF_type[i]==2)
//				//	mClr[i] = COLORA(1,1,1,1);
//			}
//			break;
//		case 1://pressure
//			for (int i=0; i<NumPoints(); i++){
//				if (mIsBound[i]==0){
//					if(mPressure[i]>100)
//						mPressure[i]=100;
//					mClr[i] = COLORA(1,1-mPressure[i]/100, 1-mPressure[i]/100, 1);
//				}
//				else {
//					mClr[i] = COLORA(1, 1, 1, 0.5);
//				}
//			}
//			break;
//		case 2://pressure force
//			for (int i=0; i<NumPoints(); i++){
//				if (mIsBound[i]==0)
//					mClr[i] = COLORA( abs(mVelEval[i].x)/10., abs(mVelEval[i].y)/10., abs(mVelEval[i].z)/10., 1);
//				else
//					mClr[i] = COLORA( 1,1,1,0.5);
//			}
//			break;
//		}
//
//
//		//Drawing   ====================================
//
//		glEnable ( GL_LIGHTING );		
//		glEnable(GL_BLEND); 
//		glEnable(GL_ALPHA_TEST); 
//		glAlphaFunc( GL_GREATER, 0.0 ); 
//		glEnable ( GL_COLOR_MATERIAL );
//		glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
//				
//		// Point sprite size
//	    glPointSize ( 64 );
//		glEnable ( GL_POINT_SIZE );		
//		
//		glEnable(GL_POINT_SPRITE_ARB); 		
//		{
//			float quadratic[] =  { 0.0f, 0.3f, 0.00f };
//			glEnable (  GL_POINT_DISTANCE_ATTENUATION  );
//			glPointParameterfvARB(  GL_POINT_DISTANCE_ATTENUATION, quadratic );
//		}
//		
//		glPointParameterfARB( GL_POINT_SIZE_MIN_ARB, 1.0f );
//
//		// Texture and blending mode
//		glEnable ( GL_TEXTURE_2D );
//		glBindTexture ( GL_TEXTURE_2D, mImg.getID() );
//		glTexEnvi (GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
//		glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_BLEND );
//		glBlendFunc ( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
//
//		// Point buffers
//		glBindBufferARB ( GL_ARRAY_BUFFER_ARB, mVBO[0] );
//		glBufferDataARB ( GL_ARRAY_BUFFER_ARB, NumPoints()*sizeof(cfloat3), mPos, GL_DYNAMIC_DRAW_ARB);		
//		glVertexPointer ( 3, GL_FLOAT, 0, 0x0 );				
//		glBindBufferARB ( GL_ARRAY_BUFFER_ARB, mVBO[1] );
//		glBufferDataARB ( GL_ARRAY_BUFFER_ARB, NumPoints()*sizeof(uint), mClr, GL_DYNAMIC_DRAW_ARB);
//		glColorPointer ( 4, GL_UNSIGNED_BYTE, 0, 0x0 ); 
//		glEnableClientState ( GL_VERTEX_ARRAY );
//		glEnableClientState ( GL_COLOR_ARRAY );
//		  
//		// Render - Point Sprites
//		glNormal3f ( 0, 1, 0.001  );
//		//glColor3f ( 1, 1, 1 );
//		glDrawArrays ( GL_POINTS, 0, NumPoints() );
//
//		// Restore state
//		glDisableClientState ( GL_VERTEX_ARRAY );
//		glDisableClientState ( GL_COLOR_ARRAY );
//		glDisable (GL_POINT_SPRITE_ARB); 
//		glDisable ( GL_ALPHA_TEST );
//		glDisable ( GL_TEXTURE_2D );
//		glDepthMask( GL_TRUE );   
//
//		break;
//	};
//}


#include <png.h>

void FluidSystem::CaptureVideo (int width, int height)
{

	
	char fileName[64];   
    sprintf( fileName, "snapshot/snapshot_%04d.png", frameNo );
 									// record frame buffer directly to image pixels
    uchar * pixels = new uchar[width*height*3];
    glReadPixels( 0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, (GLvoid*)pixels );	

    png_structp png = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
    if (!png)
        return;

    png_infop info = png_create_info_struct(png);
    if (!info){
        png_destroy_write_struct(&png, &info);
        return;
    }

    FILE* fp = fopen(fileName, "wb");
    png_init_io(png, fp);
    png_set_IHDR(png, info, width, height, 8, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_BASE, PNG_FILTER_TYPE_BASE);
    png_colorp palette = (png_colorp)png_malloc(png, PNG_MAX_PALETTE_LENGTH * sizeof(png_color));
    if (!palette){
        fclose(fp);
        png_destroy_write_struct(&png, &info);
        return;
    }
    png_set_PLTE(png, info, palette, PNG_MAX_PALETTE_LENGTH);
    png_write_info(png, info);
    png_set_packing(png);

    png_bytepp rows = (png_bytepp)png_malloc(png, height*sizeof(png_bytep));
    for (int i=0; i<height; i++)
        rows[i] = (png_bytep)(pixels + (height-i-1)*width*3);
    png_write_image(png,rows);
    png_write_end(png,info);
    png_free(png,palette);
	png_free(png,rows);
    png_destroy_write_struct(&png, &info);

    fclose(fp);
    
	delete pixels;

    return;
}

void FluidSystem::saveParticle(std::string name)
{
	TransferFromCUDAForLoad();
	FILE* fp;
	int n = 0;
	for (int i = 0;i<pointNum;i++)
		if (mIsBound[i] == 0 && mPos[i].x>-500)
			n++;
	fp = fopen(name.c_str(),"w");
	fprintf(fp,"%d\n",n);
	cfloat3* ppos = mPos;
	for (int i = 0;i<pointNum;i++,ppos++)
	if (mIsBound[i] == 0 && mPos[i].x>-500)
	{
		fprintf(fp,"%f %f %f",ppos->x,ppos->y,ppos->z);
		for (int j = 0;j<MAX_FLUIDNUM;j++)
		{
			fprintf(fp," %f",*(m_alpha+i*MAX_FLUIDNUM + j));
			fprintf(fp," %f",*(m_alpha_pre+i*MAX_FLUIDNUM + j));
		}
		fprintf(fp," %f",m_restMass[i]);
		fprintf(fp," %f",m_restDensity[i]);
		fprintf(fp," %f",m_visc[i]);
		fprintf(fp," %f %f %f",mVel[i].x, mVel[i].y, mVel[i].z);
		fprintf(fp," %f %f %f",mVelEval[i].x, mVelEval[i].y, mVelEval[i].z);
		fprintf(fp," %d",MF_type[i]);
		if( MF_type[i]==1)
			for(int k=0; k<9; k++){
				fprintf(fp," %f",MF_tensor[i*9+k]);
			}
	}
	fclose(fp);
	//SetYan(SAVE_STAT,0);
}

int FluidSystem::loadParticle(std::string name)
{
	int n,p;
	float f1,f2,f3;
	FILE* fp;
	fp = fopen(name.c_str(),"r");
	fscanf(fp,"%d",&n);
	cfloat3* ppos = mPos;
	for (int i = 0;i<n;i++)
	{
		p = AddParticle(); // the index of the added particle
		if (p!=-1)
		{
			fscanf(fp,"%f %f %f",&f1,&f2,&f3);
			(ppos+p)->Set(f1,f2,f3);
			for (int j = 0;j<MAX_FLUIDNUM;j++)
			{
				fscanf(fp,"%f",(m_alpha+p*MAX_FLUIDNUM + j));
				fscanf(fp,"%f",(m_alpha_pre+p*MAX_FLUIDNUM + j));
			}
			if(sceneID==3){
				for (int j = 0;j<MAX_FLUIDNUM;j++)
				{
					m_alpha[p*MAX_FLUIDNUM] = 0.15f;
					m_alpha[p*MAX_FLUIDNUM+1] = 0.25f;
					m_alpha[p*MAX_FLUIDNUM+2] = 0.4f;
					m_alpha[p*MAX_FLUIDNUM+3] = 0.2f;
					m_alpha_pre[p*MAX_FLUIDNUM] = 0.15f;
					m_alpha_pre[p*MAX_FLUIDNUM+1] = 0.25f;
					m_alpha_pre[p*MAX_FLUIDNUM+2] = 0.4f;
					m_alpha_pre[p*MAX_FLUIDNUM+3] = 0.2f;
				}
			}
			if(sceneID==5){
				for (int j = 0;j<MAX_FLUIDNUM;j++)
				{
					m_alpha[p*MAX_FLUIDNUM] = 0.0f;
					m_alpha[p*MAX_FLUIDNUM+1] = 0.0f;
					m_alpha[p*MAX_FLUIDNUM+2] = 0.0f;
					m_alpha[p*MAX_FLUIDNUM+3] = 1.0f;
	//				m_alpha[p*MAX_FLUIDNUM+3] = 0.2f;
	//				m_alpha[p*MAX_FLUIDNUM+4] = 0.2f;
					m_alpha_pre[p*MAX_FLUIDNUM] = 0.0f;
					m_alpha_pre[p*MAX_FLUIDNUM+1] = 0.0f;
					m_alpha_pre[p*MAX_FLUIDNUM+2] = 0.0f;
					m_alpha_pre[p*MAX_FLUIDNUM+3] = 1.0f;
	//				m_alpha_pre[p*MAX_FLUIDNUM+3] = 0.2f;
	//				m_alpha_pre[p*MAX_FLUIDNUM+4] = 0.2f;
				}
			}
			fscanf(fp,"%f",m_restMass+p);
			fscanf(fp,"%f",m_restDensity+p);
			fscanf(fp,"%f",m_visc+p);
			fscanf(fp," %f %f %f",&mVel[i].x, &mVel[i].y, &mVel[i].z);
			fscanf(fp," %f %f %f",&mVelEval[i].x, &mVelEval[i].y, &mVelEval[i].z);
			fscanf(fp,"%d",&MF_type[i]);

			if(MF_type[i]==1){
				for(int k=0; k<9; k++){
					fscanf(fp,"%f",&MF_tensor[i*9+k]);
				}
			}
			*(mClr+p) = COLORA( 1,1,1,1);
		}
	}
	fclose(fp);
	return n;
}

void FluidSystem::record ( int param, std::string name, cTime& start )
{
	//cTime stop;
	//m_Param [ param ] = stop - start;
}