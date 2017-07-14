/*
FLUIDS v.3 - SPH Fluid Simulator for CPU and GPU
Copyright (C) 2012. Rama Hoetzlein, http://www.rchoetzlein.com

ZLib license
This software is provided 'as-is', without any express or implied
warranty.  In no event will the authors be held liable for any damages
arising from the use of this software.

Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it
freely, subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not
claim that you wrote the original software. If you use this software
in a product, an acknowledgment in the product documentation would be
appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be
misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "SolverGUI.h"

bool bTiming = true;
bool bRec = false;
int mFrame = 0;

SolverGUI		solverGUI;

//Camera3D		cam;

float window_width  = 1024;
float window_height = 768;

//Vector3DF	obj_from, obj_angs, obj_dang;
//Vector4DF	light[2], light_to[2];				// Light stuff
float		light_fov;

int		psys_rate = 0;							// Particle stuff
int		psys_freq = 1;
int		psys_playback;

bool	bHelp = false;					// Toggles
int		iShade = 0;						// Shading mode (default = no shadows)
int		iClrMode = 0;
bool    bPause = false;


										// Mouse control
#define DRAG_OFF		0				// mouse states
#define DRAG_LEFT		1
#define DRAG_RIGHT		2
#define DRAG_UP			3
#define DRAG_DOWN		4
#define DRAG_MIDDLE     5
int		last_x = -1, last_y = -1;		// mouse vars
int		mode = 0;
int		dragging = 0;
int		psel;

GLuint	screen_id;
GLuint	depth_id;


// Different things we can move around
#define MODE_CAM		0
#define MODE_CAM_TO		1
#define MODE_OBJ		2
#define MODE_OBJPOS		3
#define MODE_OBJGRP		4
#define MODE_LIGHTPOS	5

#define MODE_DOF		6


GLuint screenBufferObject;
GLuint depthBufferObject;
GLuint envid;



//void SolverGUI::drawScene(float* viewmat, bool bShade)
//{
//	if (iShade <= 1 && bShade) {
//
//		glEnable(GL_LIGHTING);
//		glEnable(GL_LIGHT0);
//		glDisable(GL_COLOR_MATERIAL);
//
//		Vector4DF amb, diff, spec;
//		float shininess = 5.0;
//
//		glColor3f(1, 1, 1);
//		glLoadIdentity();
//		glLoadMatrixf(viewmat);
//
//		float pos[4];
//		pos[0] = light[0].x;
//		pos[1] = light[0].y;
//		pos[2] = light[0].z;
//		pos[3] = 1;
//		amb.Set(0, 0, 0, 1); diff.Set(1, 1, 1, 1); spec.Set(1, 1, 1, 1);
//		glLightfv(GL_LIGHT0, GL_POSITION, (float*)&pos[0]);
//		glLightfv(GL_LIGHT0, GL_AMBIENT, (float*)&amb.x);
//		glLightfv(GL_LIGHT0, GL_DIFFUSE, (float*)&diff.x);
//		glLightfv(GL_LIGHT0, GL_SPECULAR, (float*)&spec.x);
//
//		amb.Set(0, 0, 0, 1); diff.Set(.3, .3, .3, 1); spec.Set(.1, .1, .1, 1);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, (float*)&amb.x);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, (float*)&diff.x);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, (float*)&spec.x);
//		glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, (float*)&shininess);
//
//
//		//glColorMaterial ( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
//
//		glLoadMatrixf(viewmat);
//		//绘制平面
//#ifdef DRAW_PLANE
//		glBegin(GL_QUADS);
//		glNormal3f(0, 1, 0.001);
//		for (float x=-1000; x <= 1000; x += 100.0) {
//			for (float y=-1000; y <= 1000; y += 100.0) {
//				glVertex3f(x, 0.0, y);
//				glVertex3f(x+100, 0.0, y);
//				glVertex3f(x+100, 0.0, y+100);
//				glVertex3f(x, 0.0, y+100);
//			}
//		}
//		glEnd();
//#endif
//
//#ifdef DRAW_LINE
//
//#ifdef DRAW_PLANE
//		glColor3f(0.1, 0.1, 0.2);
//#else
//		glColor3f(0.6, 0.6, 0.6);
//#endif
//
//		glDisable(GL_LIGHTING);
//		glBegin(GL_LINES);
//		//绘制网格线
//		for (float n=-100; n <= 100; n += 10.0) {
//			glVertex3f(-100, 0.1, n);
//			glVertex3f(100, 0.1, n);
//			glVertex3f(n, 0.1, -100);
//			glVertex3f(n, 0.1, 100);
//		}
//		glVertex3f(light[0].x, light[0].y, 0);
//		glVertex3f(light[0].x, light[0].y, light[0].z);
//		glEnd();
//
//		psys->Draw(cam, 0.8);				// Draw particles		
//#endif
//	}
//
//}




//mint::Time nowTime;
double dst;

//void SolverGUI::draw2D()
//{
//
//	glDisable(GL_LIGHTING);
//	glDisable(GL_DEPTH_TEST);
//
//	glMatrixMode(GL_PROJECTION);
//	glLoadIdentity();
//	glScalef(2.0/window_width, -2.0/window_height, 1);		// Setup view (0,0) to (800,600)
//	glTranslatef(-window_width/2.0, -window_height/2, 0.0);
//
//
//	float view_matrix[16];
//	glMatrixMode(GL_MODELVIEW);
//	glLoadIdentity();
//
//	char disp[200];
//
//	/*psys->getModeClr ();
//	strcpy ( disp, psys->getModeStr().c_str() ); drawText ( 20, 40, disp );*/
//
//	glColor4f(1.0, 1.0, 1.0, 1.0);
//	strcpy(disp, "Press H for help.");		drawText(10, 20, disp);
//
//	glColor4f(1.0, 1.0, 0.0, 1.0);
//	strcpy(disp, "");
//	if (psys->GetToggle(PCAPTURE)) strcpy(disp, "CAPTURING VIDEO");
//	drawText(200, 20, disp);
//	if (psys->GetYan(START_OUTPUT) == 1) strcpy(disp, "OUTPUTING DATA");
//	drawText(500, 20, disp);
//	if (bHelp) {
//
//		sprintf(disp, "Mode:                %s", psys->getModeStr().c_str());					drawText(20, 40, disp);
//
//		sprintf(disp, "Scene:               %s (id: %d)", psys->getSceneName().c_str(), (int)psys->GetParam(PEXAMPLE));				drawText(20, 60, disp);
//
//		sprintf(disp, "# Particles:         %d", psys->NumPoints());					drawText(20, 80, disp);
//
//		sprintf(disp, "Grid Density:        %f", psys->GetParam(PGRID_DENSITY));		drawText(20, 100, disp);
//		
//		//sprintf(disp, "Grid Count:          %f", (float)psys->GetParam(PSTAT_GRIDCNT) / psys->GetParam(PSTAT_OCCUPY));	drawText(20, 110, disp);
//		//sprintf(disp, "Grid Occupancy:      %f%%", (float)psys->GetParam(PSTAT_OCCUPY) / psys->getGridTotal());		drawText(20, 130, disp);
//		
//		sprintf(disp, "Grid Resolution:     %d x %d x %d (%d)", (int)psys->GetGridRes().x, (int)psys->GetGridRes().y, (int)psys->GetGridRes().z, psys->getGridTotal());		drawText(20, 140, disp);
//		
//		int nsrch = pow(psys->getSearchCnt(), 1/3.0);
//		sprintf(disp, "Grid Search:         %d x %d x %d", nsrch, nsrch, nsrch);			drawText(20, 150, disp);
//		//sprintf(disp, "Search Count:        %d, ave: %f, max: %f", (int)psys->GetParam(PSTAT_SRCH), psys->GetParam(PSTAT_SRCH)/psys->NumPoints(), psys->GetParam(PSTAT_SRCHMAX)/psys->NumPoints());		drawText(20, 160, disp);
//		//sprintf(disp, "Neighbor Count:      %d, ave: %f, max: %f", (int)psys->GetParam(PSTAT_NBR), psys->GetParam(PSTAT_NBR)/psys->NumPoints(), psys->GetParam(PSTAT_NBRMAX)/psys->NumPoints());		drawText(20, 170, disp);
//		//sprintf(disp, "Search Overhead:     %.2fx", psys->GetParam(PSTAT_SRCH)/psys->GetParam(PSTAT_NBR));		drawText(20, 180, disp);
//
//		//sprintf(disp, "Insert Time:         %.3f ms", psys->GetParam(PTIME_INSERT));			drawText(20, 200, disp);
//		sprintf(disp, "Sort Time:           %.3f ms", psys->GetParam(PTIME_SORT));			drawText(20, 210, disp);
//		//sprintf(disp, "Count Time:          %.3f ms", psys->GetParam(PTIME_COUNT));			drawText(20, 220, disp);
//		sprintf(disp, "Pressure Time:       %.3f ms", psys->GetParam(PTIME_PRESS));			drawText(20, 230, disp);
//		sprintf(disp, "Force Time:          %.3f ms", psys->GetParam(PTIME_FORCE));			drawText(20, 240, disp);
//		sprintf(disp, "Advance Time:        %.3f ms", psys->GetParam(PTIME_ADVANCE));			drawText(20, 250, disp);
//
//		//multi fluid
//		float st = psys->GetParam(PTIME_SORT) + psys->GetParam(PTIME_PRESS)+psys->GetParam(PTIME_FORCE) + psys->GetParam(PTIME_ADVANCE) + psys->GetParam(PTIMEDRIFTVEL) + psys->GetParam(PTIMEALPHA) + psys->GetParam(PTIMECORR);
//
//		//for performance
//		if (!bPause) {
//			mint::Time nowTime;
//			nowTime.SetSystemTime(ACC_MSEC);
//			double now = nowTime.GetMSec();
//			if (psys->lastTime == 0) {
//				psys->lastTime = now;
//			}
//			dst = now - psys->lastTime;
//			if (dst == 0)
//				dst = 1.0f;
//			psys->lastTime = now;
//		}
//		sprintf(disp, "Time Performance:      %.3f ms, %.1f fps", dst, tm_fps);									drawText(20, 260, disp);
//		sprintf(disp, "Simulation Performance: %.3f ms", st); drawText(20, 270, disp);
//		//sprintf ( disp,	"Performance:         %d particles/sec", (int) ((psys->NumPoints()*1000.0)/dst) );			drawText ( 20, 270,  disp );
//
//		//sprintf(disp, "Particle Memory:     %.4f MB", (float)psys->GetParam(PSTAT_PMEM)/1000000.0f);		drawText(20, 290, disp);
//		//sprintf(disp, "Grid Memory:         %.4f MB", (float)psys->GetParam(PSTAT_GMEM)/1000000.0f);		drawText(20, 300, disp);
//
//		size_t fm, tm;
//		cudaMemGetInfo(&fm, &tm);
//		sprintf(disp, "GPU Free Memory:     %u MB", fm/1000000); drawText(20, 310, disp);
//		sprintf(disp, "GPU Total Memory:    %u MB", tm/1000000); drawText(20, 320, disp);
//		sprintf(disp, "GPU Memory Load:     %.2f%%",  (1-(float)fm/tm)*100); drawText(20,330, disp);
//		//multi fluid
//		
//		//sprintf(disp, "Drift Velocity Time:        %.3f ms", psys->GetParam(PTIMEDRIFTVEL));			drawText(20, 330, disp);
//		//sprintf(disp, "Alpha Advance Time:        %.3f ms", psys->GetParam(PTIMEALPHA));			drawText(20, 340, disp);
//		//sprintf(disp, "Alpha Correct Time:        %.3f ms", psys->GetParam(PTIMECORR));			drawText(20, 350, disp);
//		//sprintf(disp, "Split Time:                %.3f ms", psys->GetParam(PTIMESPLIT));			drawText(20, 360, disp);
//		sprintf(disp, "Frame NO:                  %d ", psys->frameNum());			drawText(20, 370, disp);
//	}
//}

int frame;



void SolverGUI::render(){

	if (!bPause) {
		psys->Run(window_width, window_height);

		frame++;
		//measureFPS();
	}

	glEnable(GL_DEPTH_TEST);

	glClearColor(0.4, 0.4, 0.4, 1);

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//glDisable(GL_CULL_FACE);
	//glEnable(GL_LIGHTING);

	// Draw Particles
	//drawScene(cam.getViewMatrix().GetDataF(), true);
	//draw2D();
	
	//get Timer
	camera.velMax = 100;


	//Update Camera Position
	clock_t Now = clock();
	if (LastTime == 0)
		LastTime = Now;

	camera.AdvanceCamera( (float)(Now-LastTime)/1000.f);
	ViewMatrix = camera.viewmat;
	LastTime = Now;

	//
	//camera.MoveCamera(cfloat3(0.01,0,0));
	
	
	DrawParticleSystem();
	//DrawCube();

	glutSwapBuffers();
	glutPostRedisplay();

}



void display()
{

	solverGUI.render();

}

void resize(int width, int height)
{
	// set window height and width
	window_width  = (float)width;
	window_height = (float)height;
	glViewport(0, 0, width, height);

	solverGUI.ReSize(width,height);
}

void SolverGUI::ReSize(int width, int height) {

	camera.SetProjParam(60, (float)width/height, 1.0f, 1000.0f);
	camera.ProjectionMat();

	glUseProgram(ShaderIds[0]);
	glUniformMatrix4fv(ProjectionMatrixUniformLocation, 1, GL_FALSE, camera.projmat.data);
	glUseProgram(0);
}


void keyboard_func(unsigned char key, int x, int y)
{
	solverGUI.keyDown(key);
}
void SolverGUI::keyDown(unsigned char key){

	switch (key) {
		//case 'R': case 'r': {
		//	psys->StartRecord();
		//} break;
		//case 'P': case 'p': 	
		/*switch ( psys->getMode() ) {
		case MODE_RECORD: case MODE_SIM:	psys_playback = psys->getLastRecording();	break;
		case MODE_PLAYBACK:					psys_playback--; if ( psys_playback < 0 ) psys_playback = psys->getLastRecording();	break;
		};
		psys->StartPlayback ( psys_playback );
		} break;*/
		//case 'M': case 'm': {
		//	psys->SetParam(PNUM, psys->GetParam(PNUM)*2, 4, 40000000);
		//	psys->Setup(false);
		//} break;
		//case 'N': case 'n': {
		//	psys->SetParam(PNUM, psys->GetParam(PNUM)/2, 4, 40000000);
		//	psys->Setup(false);
		//} break;
		//case '0':
		//	UpdateEmit();
		//	psys_freq++;
		//	psys->SetVec(PEMIT_RATE, Vector3DF(psys_freq, psys_rate, 0));
		//	break;
		//case '9':
		//	UpdateEmit();
		//	psys_freq--;  if (psys_freq < 0) psys_freq = 0;
		//	psys->SetVec(PEMIT_RATE, Vector3DF(psys_freq, psys_rate, 0));
		//	break;
		//case '.': case '>':
		//	UpdateEmit();
		//	if (++psys_rate > 100) psys_rate = 100;
		//	psys->SetVec(PEMIT_RATE, Vector3DF(psys_freq, psys_rate, 0));
		//	break;
		//case ',': case '<':
		//	UpdateEmit();
		//	if (--psys_rate < 0) psys_rate = 0;
		//	psys->SetVec(PEMIT_RATE, Vector3DF(psys_freq, psys_rate, 0));
		//	break;
	case 'B': case 'b':
		psys->bOutput = ! psys->bOutput;
		break;
	case '5':

		break;
	case '6':
		psys->displayMode ++;
		if(psys->displayMode > 2)
			psys->displayMode = 0;
		break;
	case '7':
		//psys->SetSimuP(SHOW_TYPE, (psys->GetSimuP(SHOW_TYPE) + 1)%2);
		break;
	case '8':
		//if (psys->GetYan(SAVE_STAT) == 0)
		//	psys->SetYan (SAVE_STAT,1);
		//psys->saveParticle("save_stat.txt");
		break;
		/*case 'f': case 'F':		psys->IncParam(PMODE, -1, 1, 8);		psys->Setup(false); break;
		case 'g': case 'G':		psys->IncParam(PMODE, 1, 1, 8);		psys->Setup(false); break;
		*/
	case ' ':				bPause = !bPause; break;//psys->SetParam ( PMODE, RUN_PAUSE, 0, 8 );	break;		// pause

	case '1':			
		//psys->IncParam(PDRAWGRID, 1, 0, 1);		
		break;
	case '2':				
		//psys->IncParam(PDRAWTEXT, 1, 0, 1);		
		break;

		/*case 'C':	mode = MODE_CAM_TO;	break;
		case 'c': 	mode = MODE_CAM;	break;*/
	case 'h': case 'H':	bHelp = !bHelp; break;
		/*case 'i': case 'I':
		UpdateEmit();
		mode = MODE_OBJPOS;
		break;
		case 'o': case 'O':
		UpdateEmit();
		mode = MODE_OBJ;
		break;*/

	case 'x': case 'X':
		if (++iClrMode > 2) iClrMode = 0;
		//psys->SetParam(PCLR_MODE, iClrMode);
		break;
	case 'l': case 'L':	mode = MODE_LIGHTPOS;	break;
	case 'j': case 'J': {
		//int d = psys->GetParam(PDRAWMODE) + 1;
		//if (d > 2) d = 0;
		//psys->SetParam(PDRAWMODE, d);
	} break;
		//case 'k': case 'K':	if (++iShade > 3) iShade = 0;		break;

	case 'a': case 'A':		
	case 'd': case 'D':		
	case 'w': case 'W':		
	case 's': case 'S':		
	case 'q': case 'Q':		
	case 'z': case 'Z':		
		camera.SetVelocty(key);
		break;

	case 27:			    exit(0); break;

	case '`':				
		psys->bSnapshot = ! psys->bSnapshot;
		break;

		//case 't': case 'T':		psys->Setup(true); break;

		/*case '-':  case '_':
		psys->IncParam(PGRID_DENSITY, -0.2, 1, 10);
		psys->Setup(true);
		break;
		case '+': case '=':
		psys->IncParam(PGRID_DENSITY, 0.2, 1, 10);
		psys->Setup(true);
		break;*/
		/*case '[':
		psys->IncParam(PEXAMPLE, -1, 0, 10);
		psys->Setup(true);
		UpdateEmit();
		break;
		case ']':
		psys->IncParam(PEXAMPLE, +1, 0, 10);
		psys->Setup(true);
		UpdateEmit();
		break;*/
	default:
		break;
	}
}



void keyUpFunc(unsigned char key, int x,int y){
	solverGUI.keyUp(key);	
}
void SolverGUI::keyUp(unsigned char key){
	switch (key) {
	case 'a': case 'A':
	case 'd': case 'D':
	case 'w': case 'W':
	case 's': case 'S':
	case 'q': case 'Q':
	case 'z': case 'Z':
		camera.UnSetVelocity(key);
		break;
	}
}


//Vector3DF cangs;
//Vector3DF ctp;
float cdist;

#define GLUT_WHEEL_UP	3      
#define GLUT_WHEEL_DOWN 4  

void mouse_click_func(int button, int state, int x, int y)
{
	if (state == GLUT_DOWN) {
		
		if (button == GLUT_LEFT_BUTTON)		
			dragging = DRAG_LEFT;
		else if (button == GLUT_RIGHT_BUTTON) 
			dragging = DRAG_RIGHT;
		else if (button == GLUT_MIDDLE_BUTTON)
			dragging = DRAG_MIDDLE;


		else if(button == GLUT_WHEEL_UP){
			
		}
		else if(button == GLUT_WHEEL_DOWN){
			
		}

	}
	else if (state==GLUT_UP) { //wheel up and down
		if (button == GLUT_WHEEL_UP)
			dragging = DRAG_UP;
		else if (button == GLUT_WHEEL_DOWN)
			dragging = DRAG_DOWN;
		else{
			dragging = DRAG_OFF;
		}
	}

	last_x = x;
	last_y = y;

	solverGUI.MouseClick(x,y,state,button);
}


void mouse_drag_func(int x,int y){
	int dx = x - last_x;
	int dy = y - last_y;
	last_x = x;
	last_y = y;
	solverGUI.drag(dx,dy);
}

void SolverGUI::MouseClick(int x,int y,int state,int button){
	if(state==GLUT_DOWN){
		if(button == GLUT_WHEEL_UP || button == GLUT_WHEEL_DOWN){
			if(button == GLUT_WHEEL_UP)
				camera.zoomin = camera.zoomin * 0.8 + 0.2 * 50;
			else
				camera.zoomin = camera.zoomin * 0.8 - 0.2 * 50;
			camera.fovy += camera.zoomin * 0.5;
			camera.fovy = min(100, max(10, camera.fovy));
			camera.ProjectionMat();
		}
	}
}

void SolverGUI::drag(int dx,int dy){

	//small displacement increment scheme
	if(dragging == DRAG_LEFT){
		float alpha = 0.4;
		camera.rotatexy = camera.rotatexy *(1-alpha) + cfloat2(dx,dy)*alpha;
	}else if(dragging == DRAG_RIGHT){
		camera.zoomin = camera.zoomin * 0.8 + 0.2 * dy;
		camera.fovy += camera.zoomin * 0.5;
		camera.fovy = min(100, max(10, camera.fovy));
		camera.ProjectionMat();
	}else if(dragging == DRAG_MIDDLE){
		float alpha = 0.4;
		camera.transxy = camera.transxy *(1-alpha) + cfloat2(dx, dy)*alpha;
	}
}


void idle_func()
{
	glutPostRedisplay();
}

//void initParameters()
//{
//	srand(time(0x0));
//	
//	// Initialize camera
//	cam.setOrbit(Vector3DF(200, 30, 0), Vector3DF(2, 30, 2), 400, 400);
//	cam.setFov(35);
//	cam.updateMatricies();
//
//	light[0].x = 0;		light[0].y = 200;	light[0].z = 0; light[0].w = 1;
//	light_to[0].x = 0;	light_to[0].y = 0;	light_to[0].z = 0; light_to[0].w = 1;
//
//	light[1].x = 55;		light[1].y = 140;	light[1].z = 50;	light[1].w = 1;
//	light_to[1].x = 0;	light_to[1].y = 0;	light_to[1].z = 0;		light_to[1].w = 1;
//
//	light_fov = 45;
//}


#include "GL/glew.h"
#include "GL/freeglut.h"

void Destroy() {
	solverGUI.DestroyCube();
}

void SolverGUI::Initialize(int argc, char** argv){
	// Initialize CUDA
	cudaInit(argc, argv);

	// set up the window
	glutInit(&argc, &argv[0]);

	glutInitContextVersion(4,0);
	glutInitContextFlags(GLUT_FORWARD_COMPATIBLE);
	glutInitContextProfile(GLUT_CORE_PROFILE);

	glutSetOption(GLUT_ACTION_ON_WINDOW_CLOSE, GLUT_ACTION_GLUTMAINLOOP_RETURNS);

	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH |GLUT_MULTISAMPLE);
	//glutInitWindowPosition(100, 100);
	glutInitWindowSize((int)window_width, (int)window_height);
	
	glutCreateWindow("Multiphase Nekobus");

	// callbacks
	glutDisplayFunc(display);
	glutReshapeFunc(resize);
	glutCloseFunc(Destroy);

	// keyboard
	glutSetKeyRepeat(GLUT_KEY_REPEAT_OFF);
	glutKeyboardFunc(keyboard_func);
	glutKeyboardUpFunc(keyUpFunc);

	glutMouseFunc(mouse_click_func);
	glutMotionFunc(mouse_drag_func);
	//glutPassiveMotionFunc(mouse_move_func);
	glutIdleFunc(idle_func);

	GLenum GlewInitResult = glewInit();
	if(GLEW_OK != GlewInitResult){
		fprintf(stderr, "ERROR: %s\n", glewGetErrorString(GlewInitResult));
		exit(EXIT_FAILURE);
	}

	fprintf(stdout, "INFO: OpenGL Version: %s\n", glGetString(GL_VERSION));
	glClearColor(0,0,0,0);

	//Camera Parameters
	ModelMatrix = IDENTITY_MAT;
	ProjectionMatrix = IDENTITY_MAT;
	camera.forceup = cfloat3(0, 1, 0);
	camera.lookat(cfloat3(0, 40, 60), cfloat3(0, 20, -40));
	ViewMatrix = camera.viewmat;
	
	//Shaders & VAO, VBO
	CreateParticleSystem();

	//initParameters();

	//Texture
	int width,height;
	unsigned char* image = SOIL_load_image("drop.png",&width,&height,0,SOIL_LOAD_RGBA);
	GLuint texture;
	glGenTextures(1, &texture);
	glBindTexture(GL_TEXTURE_2D,texture);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, image);
	glGenerateMipmap(GL_TEXTURE_2D);
	SOIL_free_image_data(image);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

}


void SolverGUI::CreateCube(){
	const cvertex VERTICES[8] =
	{
		{{-.5f, -.5f,  .5f, 1},{0, 0, 1, 1}},
		{{-.5f,  .5f,  .5f, 1},{1, 0, 0, 1}},
		{{.5f,  .5f,  .5f, 1},{0, 1, 0, 1}},
		{{.5f, -.5f,  .5f, 1},{1, 1, 0, 1}},
		{{-.5f, -.5f, -.5f, 1},{1, 1, 1, 1}},
		{{-.5f,  .5f, -.5f, 1},{1, 0, 0, 1}},
		{{.5f,  .5f, -.5f, 1},{1, 0, 1, 1}},
		{{.5f, -.5f, -.5f, 1},{0, 0, 1, 1}}
	};
	const GLuint INDICES[36] =
	{
		0,2,1,  0,3,2,
		4,3,0,  4,7,3,
		4,1,5,  4,0,1,
		3,6,2,  3,7,6,
		1,6,5,  1,2,6,
		7,5,6,  7,4,5
	};

	ShaderIds[0] = glCreateProgram();
	ExitOnGLError("ERROR: Could not create the shader program");

	ShaderIds[1] = LoadShader("shader/SimpleShader.fragment.glsl", GL_FRAGMENT_SHADER);
	ShaderIds[2] = LoadShader("shader/SimpleShader.vertex.glsl", GL_VERTEX_SHADER);
	ShaderIds[3] = LoadShader("shader/SimpleShader.geometry.glsl", GL_GEOMETRY_SHADER);
	glAttachShader(ShaderIds[0], ShaderIds[1]);
	glAttachShader(ShaderIds[0], ShaderIds[2]);
	glAttachShader(ShaderIds[0], ShaderIds[3]);

	glLinkProgram(ShaderIds[0]);
	ExitOnGLError("ERROR: Could not link the shader program");

	ModelMatrixUniformLocation = glGetUniformLocation(ShaderIds[0], "ModelMatrix");
	ViewMatrixUniformLocation = glGetUniformLocation(ShaderIds[0], "ViewMatrix");
	ProjectionMatrixUniformLocation = glGetUniformLocation(ShaderIds[0], "ProjectionMatrix");
	ParticleSizeUniformLocation = glGetUniformLocation(ShaderIds[0],"particle_size");
	ExitOnGLError("ERROR: Could not get the shader uniform locations");

	glGenBuffers(2, &BufferIds[1]);
	ExitOnGLError("ERROR: Could not generate the buffer objects");

	glGenVertexArrays(1, &BufferIds[0]);
	ExitOnGLError("ERROR: Could not generate the VAO");
	glBindVertexArray(BufferIds[0]);
	ExitOnGLError("ERROR: Could not bind the VAO");

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	ExitOnGLError("ERROR: Could not enable vertex attributes");

	glBindBuffer(GL_ARRAY_BUFFER, BufferIds[1]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(VERTICES), VERTICES, GL_STATIC_DRAW);
	ExitOnGLError("ERROR: Could not bind the VBO to the VAO");

	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, sizeof(VERTICES[0]), (GLvoid*)0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(VERTICES[0]), (GLvoid*)sizeof(VERTICES[0].position));
	ExitOnGLError("ERROR: Could not set VAO attributes");

	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BufferIds[2]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(INDICES), INDICES, GL_STATIC_DRAW);
	ExitOnGLError("ERROR: Could not bind the IBO to the VAO");

	glBindVertexArray(0);


	//Rendering Cubes
	//glEnable(GL_DEPTH_TEST);
	//glDepthFunc(GL_LESS);
	//glEnable(GL_CULL_FACE);
	//glCullFace(GL_BACK);
	//glFrontFace(GL_CCW);
	
}

void SolverGUI::CreateParticleSystem() {

	float pos[] = {-.5,-.5,.5,
		-.5,.5,.5,
		.5,.5,.5,
		.5,-.5,.5,
		-.5,-.5,-.5,
		-.5,.5,-.5,
		.5,.5,-.5,
		.5,-.5,-.5};
	float color[]={ 0,0,1,1,
		1,0,0,1,
		0,1,0,1,
		1,1,0,1,
		1,1,1,1,
		1,0,0,1,
		1,0,1,1,
		0,0,1,1};
	displayPack tmpbuf[]={
		{ cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)},
		{cfloat3(-.5,-.5,.5) , cfloat4(0,0,1,1)}
	};

	ShaderIds[0] = glCreateProgram();
	ExitOnGLError("ERROR: Could not create the shader program");

	ShaderIds[1] = LoadShader("shader/SimpleShader.fragment.glsl", GL_FRAGMENT_SHADER);
	ShaderIds[2] = LoadShader("shader/SimpleShader.vertex.glsl", GL_VERTEX_SHADER);
	ShaderIds[3] = LoadShader("shader/SimpleShader.geometry.glsl", GL_GEOMETRY_SHADER);
	glAttachShader(ShaderIds[0], ShaderIds[1]);
	glAttachShader(ShaderIds[0], ShaderIds[2]);
	glAttachShader(ShaderIds[0], ShaderIds[3]);

	glLinkProgram(ShaderIds[0]);
	ExitOnGLError("ERROR: Could not link the shader program");

	ModelMatrixUniformLocation = glGetUniformLocation(ShaderIds[0], "ModelMatrix");
	ViewMatrixUniformLocation = glGetUniformLocation(ShaderIds[0], "ViewMatrix");
	ProjectionMatrixUniformLocation = glGetUniformLocation(ShaderIds[0], "ProjectionMatrix");
	ParticleSizeUniformLocation = glGetUniformLocation(ShaderIds[0], "particle_size");
	ExitOnGLError("ERROR: Could not get the shader uniform locations");

	glGenBuffers(2, &BufferIds[1]);
	ExitOnGLError("ERROR: Could not generate the buffer objects");

	glGenVertexArrays(1, &BufferIds[0]);
	ExitOnGLError("ERROR: Could not generate the VAO");
	glBindVertexArray(BufferIds[0]);
	ExitOnGLError("ERROR: Could not bind the VAO");

	glEnableVertexAttribArray(0);
	glEnableVertexAttribArray(1);
	ExitOnGLError("ERROR: Could not enable vertex attributes");

	maxPointNum = 20000; //100,000 needed for SPH
	pointNum = 8;
	glBindBuffer(GL_ARRAY_BUFFER, BufferIds[1]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(cfloat3)*maxPointNum, NULL, GL_STATIC_DRAW);
	//glBufferSubData(GL_ARRAY_BUFFER,0,sizeof(cfloat3)*pointNum,pos);
	//glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, (GLvoid*)0);
	glBufferData(GL_ARRAY_BUFFER, sizeof(displayPack)*maxPointNum, NULL, GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0 ,sizeof(displayPack)*pointNum, tmpbuf);
	glVertexAttribPointer(0,3,GL_FLOAT,GL_FALSE,sizeof(displayPack),0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, sizeof(displayPack), (void*)sizeof(cfloat3));

	//glBindBuffer(GL_ARRAY_BUFFER, BufferIds[2]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*4*maxPointNum, NULL, GL_STATIC_DRAW);
	//glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(color), color);
	//glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, 0, 0);
	
	//glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, BufferIds[2]);
	//glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(INDICES), INDICES, GL_STATIC_DRAW);
	//ExitOnGLError("ERROR: Could not bind the IBO to the VAO");

	glBindVertexArray(0);

	
}

void SolverGUI::DestroyCube(){
	glDetachShader(ShaderIds[0], ShaderIds[1]);
	glDetachShader(ShaderIds[0], ShaderIds[2]);
	glDeleteShader(ShaderIds[1]);
	glDeleteShader(ShaderIds[2]);
	glDeleteProgram(ShaderIds[0]);
	ExitOnGLError("ERROR: Could not destroy the shaders");

	glDeleteBuffers(2, &BufferIds[1]);
	glDeleteVertexArrays(1, &BufferIds[0]);
	ExitOnGLError("ERROR: Could not destroy the buffer objects");
}

float pos[] ={-.5,-.5,.5,
-.5,.5,.5,
.5,.5,.5,
.5,-.5,.5,
-.5,-.5,-.5,
-.5,.5,-.5,
.5,.5,-.5,
.5,-.5,-.5};
float color[]={0,0,1,1,
1,0,0,1,
0,1,0,1,
1,1,0,1,
1,1,1,1,
1,0,0,1,
1,0,1,1,
0,0,1,1};

void SolverGUI::DrawParticleSystem(){
	
	

	ModelMatrix = IDENTITY_MAT;

	glUseProgram(ShaderIds[0]);
	ExitOnGLError("ERROR: Could not use the shader program");

	glUniformMatrix4fv(ModelMatrixUniformLocation, 1, GL_FALSE, ModelMatrix.data);
	glUniformMatrix4fv(ViewMatrixUniformLocation, 1, GL_FALSE, ViewMatrix.data);
	
	glUniformMatrix4fv(ProjectionMatrixUniformLocation, 1, GL_FALSE, camera.projmat.data);
	
	glUniform1f(ParticleSizeUniformLocation, 1.0f);
	ExitOnGLError("ERROR: Could not set the shader uniforms");

	glBindVertexArray(BufferIds[0]);
	ExitOnGLError("ERROR: Could not bind the VAO for drawing purposes");

	pointNum = psys->pointNum;
	//Update Particle Data
	glBindBuffer(GL_ARRAY_BUFFER, BufferIds[1]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*3*maxPointNum, NULL, GL_STATIC_DRAW);
	//glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*3*pointNum, psys->mPos);
	glBufferData(GL_ARRAY_BUFFER, sizeof(displayPack)*maxPointNum, NULL,GL_STATIC_DRAW);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(displayPack)*pointNum, psys->displayBuffer);

	//glBindBuffer(GL_ARRAY_BUFFER, BufferIds[2]);
	//glBufferData(GL_ARRAY_BUFFER, sizeof(float)*4*maxPointNum, NULL, GL_STATIC_DRAW);
	//glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float)*4*pointNum, psys->mColor);
	
	glDrawArrays(GL_POINTS, 0, pointNum);
	ExitOnGLError("ERROR: Could not draw the cube");

	glBindVertexArray(0);
	glUseProgram(0);
}

void SolverGUI::DrawCube(){
	float CubeAngle;
	clock_t Now = clock();
	if (LastTime == 0)
		LastTime = Now;

	CubeRotation += 45.0f * ((float)(Now - LastTime) / CLOCKS_PER_SEC);
	CubeAngle = deg2rad(CubeRotation);
	LastTime = Now;

	ModelMatrix = IDENTITY_MAT;
	RotateAboutY(ModelMatrix, CubeAngle);
	RotateAboutX(ModelMatrix, CubeAngle);

	glUseProgram(ShaderIds[0]);
	ExitOnGLError("ERROR: Could not use the shader program");

	glUniformMatrix4fv(ModelMatrixUniformLocation, 1, GL_FALSE, ModelMatrix.data);
	glUniformMatrix4fv(ViewMatrixUniformLocation, 1, GL_FALSE, ViewMatrix.data);
	glUniform1f(ParticleSizeUniformLocation, 0.05f);
	ExitOnGLError("ERROR: Could not set the shader uniforms");

	glBindVertexArray(BufferIds[0]);
	ExitOnGLError("ERROR: Could not bind the VAO for drawing purposes");

	//glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_INT, (GLvoid*)0);
	//glPointSize(10);
	//glBindTexture(GL_TEXTURE_2D, img.getID());

	glDrawArrays(GL_POINTS, 0, 8);
	ExitOnGLError("ERROR: Could not draw the cube");

	glBindVertexArray(0);
	glUseProgram(0);
}

void SolverGUI::SetupSolver(){
	psys = new FluidSystem();
	psys->Setup(true);
	//psys->SetupRender();
}

void SolverGUI::Run(){
	glutMainLoop();
}

void SolverGUI::Exit(){
	psys->Exit();
	delete psys;
}