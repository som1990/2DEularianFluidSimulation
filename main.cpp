/*Author: Soumitra Goswami
//Date: 4/8/2018
//Description: This is a test on a real time CPU Based 2D Eularian Fluid simulation. I test out different advection 
//             schemes including Semi Lagrangian, BFECC and Modified MacCormack. 
//             I utilize the Gauss-Seidel Pressure projection to solve incompressability.
//             Here you can paint the density to begin the simulation. The program comes with a source advection at bottom 
//             of the screen.   

*/


#include <Windows.h>
#include <omp.h>
#include <GL\glew.h>
#include <GL\freeglut.h>
#include <glm/glm.hpp>
#include <iostream>
#include <string>

#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image.h"
#include "stb_image_write.h"
#include "FluidSim.h"

using namespace std;

int iWidth, iHeight, iChannels;
float *pixmap;
float *imageFile;
unsigned char *outputFile;

int frame;
float brightnessScale;

int paint_mode, display_mode;
enum { PAINT_OBSTRUCTION, PAINT_SOURCE, PAINT_DIVERGENCE, PAINT_COLOR };
enum { COLOR_DISPLAY, DENSITY_DISPLAY, VELOCITY_DISPLAY, PRESSURE_DISPLAY, DIVERGENCE_DISPLAY };

#define SWAP(x,y) {float *temp=x;x=y;y=temp; }
#define BRUSH_SIZE 21
float source_brush[BRUSH_SIZE][BRUSH_SIZE];
float obstruction_brush[BRUSH_SIZE][BRUSH_SIZE];
int nXGrids, nYGrids;
float xLength, yLength;
float dXSize, dYSize;

int x_mouse_prev, y_mouse_prev;

float timeStep = 0.05f;

float *source_density, *source_obstruction;

//INPUTS
float gravity[2] = { 0,40.0f };
float gauss_sidel_Steps = 35;
float iterative_Orthogonal_ProjectionSteps = 1;

float kinematic_viscocity = 0;//.0001;// Amount of viscocity.
float ST_Strength = 0;// glm::exp(-0.5);// glm::exp(-1);//Amount of surface Tension
float vorticity_strength = 500;


FluidSim sim(gravity, kinematic_viscocity, ST_Strength, vorticity_strength, gauss_sidel_Steps, iterative_Orthogonal_ProjectionSteps, FluidSim::Advection::MM);

bool pause_sim = false;
bool capture_screen = false;
//-------------------------------IO FUNCTIONS------------------------------------//

void readImage(const char* fileName, float *&map)
{
	int xres, yres, channels;
	float *inImage = stbi_loadf(fileName, &xres, &yres, &channels, 0);
	if (!inImage) { return; }
	iWidth = xres;
	iHeight = yres;
	iChannels = channels;
	map = new float[xres*yres*channels];
	long index = 0;
	for (int j = 0; j < yres; j++)
	{
		for (int i = 0; i < xres; i++)
		{
			for (int c = 0; c < channels; c++)
			{
				long iIndex = (i + xres*(yres - j - 1))*channels + c;
				//cout << inImage[index] << endl;
				(map)[iIndex] = inImage[index++];
				//cout << map[iIndex] << endl;
			}
		}
	}
	//delete[] inImage;
}

void writeImage(const char* fileName, float *&map)
{
	int quality = 100;
	for (int j = 0; j < iHeight; j++)
	{
		for (int i = 0; i < iWidth; i++)
		{
			int index = i + j*(iWidth);
			int out_index = (iHeight - j - 1) * iWidth + i;
			int r = int(map[index * 3 + 0] * 255);
			int g = int(map[index * 3 + 1] * 255);
			int b = int(map[index * 3 + 2] * 255);
			
			if (r < 0) r = 0; if (r > 255) r = 255;
			if (g < 0) g = 0; if (g > 255) g = 255;
			if (b < 0) b = 0; if (b > 255) b = 255;
			
			outputFile[out_index * 3 + 0] = r;
			outputFile[out_index * 3 + 1] = g;
			outputFile[out_index * 3 + 2] = b;
		}
	}

	stbi_write_jpg(fileName, iWidth, iHeight, 3, outputFile, quality);
	

}

//--------------------------------------------INITIALIZE FUNCTIONS---------------------------------------------------//

void initMaps(float(*&map), int size, float value)
{

#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		(map)[i] = value;
	}
}

void initializeBrush()
{
	int brush_rad = (BRUSH_SIZE - 1) / 2;
	for (int j = -brush_rad; j <= brush_rad; j++)
	{
		int jj = j + brush_rad;
		float jRatio = (float(brush_rad) - fabs(j)) / float(brush_rad);
		for (int i = -brush_rad; i <= brush_rad; i++)
		{
			int ii = i + brush_rad;
			float iRatio = (float(brush_rad) - fabs(i)) / float(brush_rad);
			float radius2 = (jRatio*jRatio + iRatio*iRatio) / 2.0;
			source_brush[ii][jj] = pow(radius2, 0.5);
			obstruction_brush[ii][jj] = 1.0 - pow(radius2, 1.0 / 4.0);
		}
	}
}

void initializeFields()
{
	int size = iWidth*iHeight;
	initMaps(source_density, size, 0.0);
	initMaps(source_obstruction, size, 1.0);

}

void printMap(float *&map, int x, int y, int c, char* name)
{
	for (int j = 0; j < y; j++)
	{
		for (int i = 0; i < x; i++)
		{
			for (int k = 0; k < c; k++)
			{
				long index = (i + x*j)*c + k;
				cout << name << " Pixel(" << i << "," << j << "): " << map[index] << endl;
			}
		}
	}
}

void setNbCores(int nb)
{
	omp_set_num_threads(nb);
}

// -------------------------------------IMAGE AND PAINT FUNCTIONS------------------------------------//

void displayImage(void)
{
	if (frame != 1)
	{
		//	printMap(old_r, iWidth, iHeight, 1);
	}
	for (int j = 0; j < iHeight; j++)
	{
		#pragma omp parallel for
		for (int i = 0; i < iWidth; i++)
		{
			int index = (i + iWidth*j) * 3;
			float r, g, b;
			r = g = b = 0;
			if (display_mode == COLOR_DISPLAY)
			{
				r = sim.resImage[index];
				g = sim.resImage[index + 1];
				b = sim.resImage[index + 2];
			}
			else if (display_mode == DENSITY_DISPLAY)
			{
				r = sim.densityField[index / 3];
				g = sim.densityField[index / 3];
				b = sim.densityField[index / 3];
			}
			else if (display_mode == VELOCITY_DISPLAY)
			{
				int i = (index / 3) * 2;
				r = fabs(sim.velField[i]);
				g = fabs(sim.velField[i + 1]);
				b = 0.0;
			}
			else if (display_mode == PRESSURE_DISPLAY)
			{
				r = 0;
				g = sim.pressureField[index / 3];
				b = 0;
			}
			else if (display_mode == DIVERGENCE_DISPLAY)
			{
				r = 0;
				g = 0;
				b = sim.divergenceField[index / 3] * 100;
			}
			pixmap[index] = r * brightnessScale;
			pixmap[index + 1] = g * brightnessScale;
			pixmap[index + 2] = b * brightnessScale;
			//cout << "Pixel: " << (index / 3) << " r: " << r << " g: " << g << " b: " << b << endl;
		}
	}
}

void scaleBrightness(float value)
{
	brightnessScale *= value;
	cout << "BRIGHTNESS: " << brightnessScale << endl;
}

void paintScreen(int x, int y)
{
	int brush_radius = (BRUSH_SIZE - 1) / 2;
	int xstart = x - brush_radius;
	int ystart = y - brush_radius;

	if (xstart < 0) { xstart = 0; }
	if (ystart < 0) { ystart = 0; }

	int xend = x + brush_radius;
	int yend = y + brush_radius;
	if (xend >= iWidth) { xend = iWidth - 1; }
	if (yend >= iHeight) { yend = iHeight - 1; }

	if (paint_mode == PAINT_SOURCE)
	{
		for (int ix = xstart; ix <= xend; ix++)
		{
			for (int iy = ystart; iy <= yend; iy++)
			{
				int index = ix + iWidth*(iHeight - iy - 1);
				imageFile[3 * index + 0] += source_brush[ix - xstart][iy - ystart] * 0.1;
				imageFile[3 * index + 1] += source_brush[ix - xstart][iy - ystart] * 0.1;
				imageFile[3 * index + 2] += source_brush[ix - xstart][iy - ystart] * 0.1;
				source_density[index] += source_brush[ix - xstart][iy - ystart];

			}
		}
	}
	if (paint_mode == PAINT_OBSTRUCTION)
	{
		for (int ix = xstart; ix <= xend; ix++)
		{
			for (int iy = ystart; iy <= yend; iy++)
			{
				int index = ix + iWidth*(iHeight - iy - 1);
				imageFile[3 * index + 0] *= obstruction_brush[ix - xstart][iy - ystart];
				imageFile[3 * index + 1] *= obstruction_brush[ix - xstart][iy - ystart];
				imageFile[3 * index + 2] *= obstruction_brush[ix - xstart][iy - ystart];
				source_obstruction[index] *= glm::clamp(obstruction_brush[ix - xstart][iy - ystart], 0.0f, 1.0f);
			}
		}
	}
	return;
}





// ---------------------------GLUT FUNCTIONS ---------------------------------------------------------//

void gRender(void)
{
	glClear(GL_COLOR_BUFFER_BIT);
	glDrawPixels(iWidth, iHeight, GL_RGB, GL_FLOAT, pixmap);
	glutSwapBuffers();
}


void gIdleState(void)
{
	displayImage();
	
	if (!pause_sim) {
		paintScreen(256, 400); // you can comment this line out to remove the source. Or change the position to move the source
		sim.simulate(imageFile, source_density, source_obstruction, timeStep);
		imageFile = sim.resImage;
		initMaps(source_density, (iWidth*iHeight), 0);
	}
	glutPostRedisplay();
	if (capture_screen)
	{
		string advection;
		if (sim.advectType == FluidSim::Advection::SL) advection = "SL";
		else if (sim.advectType == FluidSim::Advection::BFECC) advection = "BFECC";
		else if (sim.advectType == FluidSim::Advection::MM) advection = "MM";
		string dispframe = to_string(frame);
		if (frame < 1000) { dispframe = "0" + dispframe; }
		if (frame < 100)  { dispframe = "0" + dispframe; }
		if (frame < 10)   { dispframe = "0" + dispframe; }
		string fName = "SG_Method_"  + advection + "_" + dispframe + ".jpg";
		
		writeImage(fName.c_str(), pixmap);
		cout << "Writing Frame: " << dispframe << endl;
	}

	frame++;

}


void gKeyboardControls(unsigned char key, int x, int y)

{
	switch (key)
	{
	case '-': case '_':
		scaleBrightness(0.9f);
		break;

	case '+': case '=':
		scaleBrightness(1.0f / 0.9f);
		break;
	case '[': case '{':
		timeStep += 0.05f;
		cout << "timeStep: " << timeStep << endl;
		break;
	case ']': case '}':
		timeStep -= 0.05f;
		cout << "timeStep: " << timeStep << endl;
		break;
	case ',': case '<':
		gravity[1] *= 0.5;
		cout << "gravity: " << gravity[1] << endl;
		break;
	case '.': case '>':
		gravity[1] *= 2;
		cout << "gravity: " << gravity[1] << endl;
		break;
	case 'c': case 'C':
		if (capture_screen == true) {
			capture_screen = false;
			cout << "CAPTURE OFF" << endl;
		}
		else {
			capture_screen = true;
			cout << "CAPTURE ON" << endl;
		}
		break;
	case 'r':
		brightnessScale = 1.0;
		break;
	case 'm': case 'M':
		display_mode++;
		display_mode = display_mode % 5;
		if		(display_mode == COLOR_DISPLAY) cout << "COLOR MODE SELECTED" << endl;
		else if (display_mode == DENSITY_DISPLAY) cout << "DENSITY MODE SELECTED" << endl;
		else if (display_mode == VELOCITY_DISPLAY) cout << "VELOCITY MODE SELECTED" << endl;
		else if (display_mode == PRESSURE_DISPLAY) cout << "PRESSURE MODE SELECTED" << endl;
		else if (display_mode == DIVERGENCE_DISPLAY) cout << "DIVERGENCE MODE SELECTED" << endl;
		break;
	case ' ':
		if (pause_sim)
		{
			cout << "SIMULATION STARTED" << endl;
			pause_sim = false;
		}
		else
		{
			cout << "SIMULATION PAUSED" << endl;
			pause_sim = true;
		}
		break;
	case 's':
		paint_mode = PAINT_SOURCE;
		cout << "PAINTING SOURCE" << endl;
		break;
	case 'o': case 'O':
		paint_mode = PAINT_OBSTRUCTION;
		cout << "PAINTING OBSTRUCTION" << endl;
		break;

	default:
		break;

	}

}

void gMouseDown(int button, int state, int x, int y)
{
	if (button != GLUT_LEFT_BUTTON) { return; }
	if (state != GLUT_DOWN) { return; }
	x_mouse_prev = x;
	y_mouse_prev = y;
	paintScreen(x, y);
}

void gMouseMove(int x, int y)
{
	x_mouse_prev = x;
	y_mouse_prev = y;
	paintScreen(x, y);
}

void PrintUsage()
{
	cout << "FLUID_paint keyboard choices\n";
	cout << "s			turns on painting source strength\n";
	cout << "o			turns on painting obstructions\n";
	cout << "+/-		increase/decrease brightness of display\n";
	cout << "r			resets brightness to default\n";
	cout << "c			toggles screen capture on/off\n";
	cout << "'['or ']'  increase/decrease timestep\n";
	cout << "',' or '.' increase/decrease Gravity\n";
	cout << "'m' or 'M' change debug mode(Color, Density, Velocity, Pressure)\n";
}

int glutFunctions(int argc, char *argv[])
{
	// Initialize GLUT
	glutInit(&argc, argv);
	// Set up some memory buffers for our display
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGBA);

	// Set the window size
	glutInitWindowSize(iWidth, iHeight);
	// Create the window with the title "SG Fluid"
	glutCreateWindow("SG Fluid");

	glClearColor(1, 1, 1, 1);
	glutDisplayFunc(&gRender);
	glutIdleFunc(&gIdleState);
	glutKeyboardFunc(&gKeyboardControls);
	glutMouseFunc(&gMouseDown);
	glutMotionFunc(&gMouseMove);

	GLenum err = glewInit();
	if (GLEW_OK != err) {
		fprintf(stderr, "GLEW error");
		return 1;
	}

	glutMainLoop();
	return 0;
}

int main(int argc, char* argv[]) {
	frame = 1;
	brightnessScale = 1;
	setNbCores(4);
	readImage("grumpy.jpg", imageFile);

	nXGrids = iWidth - 2;
	nYGrids = iHeight - 2;

	xLength = (float)iWidth;
	yLength = (float)iHeight;


	dXSize = xLength / float(nXGrids + 2);
	dYSize = yLength / float(nYGrids + 2);


	int size = iWidth*iHeight;
	int grids = (nXGrids + 2)*(nYGrids + 2);

	pixmap = new float[size * 3];
	cout << "Initializing pixmap" << endl;
	outputFile = new unsigned char[size * 3];

	//initMaps(&pixmap, size * 3, 0.0);
	source_density = new float[size];
	source_obstruction = new float[size];
	initializeFields();
	sim.initFields(nXGrids, nYGrids, dXSize, dYSize);



	initializeBrush();
	paint_mode = PAINT_SOURCE;
	display_mode = COLOR_DISPLAY;
	PrintUsage();

	int r = glutFunctions(argc, argv);

	return r;
}