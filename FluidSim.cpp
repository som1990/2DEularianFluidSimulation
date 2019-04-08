#include "FluidSim.h"



void initMaps(float (*&map), int size, float value, size_t vecSize)
{
	#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		int index = 0;
		if (vecSize == 1)
		{
			index = i;
			(map)[i] = value;
		}

		else if (vecSize == 2)
		{
			index = i * 2;
			(map)[i    ] = value;
			(map)[i + 1] = value;
		}
		else if (vecSize == 3)
		{
			index = i * 3;
			(map)[i    ] = value;
			(map)[i + 1] = value;
			(map)[i + 2] = value;
 		}

	}
	
}


void FluidSim::initFields(int nGridsX, int nGridsY, float dXSize, float dYSize)
{
	this->nGridX = nGridsX;
	this->nGridY = nGridsY;
	this->dXSize = dXSize;
	this->dYSize = dYSize;

	int size = (nGridsX + 2)*(nGridsY + 2);

	//1D Float fields
	densityField = new float[size];
	pressureField = new float[size];
	divergenceField = new float[size];
	vorticityField = new float[size];

	//2D Velocity Field
	velField = new float[size * 2];

	//2D Characterictic Maps
	characteristicMap = new float[size * 2];
	BFECC_charMap = new float[size * 2];
	MM_charMap = new float[size * 2];

	//3D Color Fields
	colorField = new float[size * 3];
	resImage = new float[size * 3];

	//Initialize Fields
	initMaps(densityField, size, 0.0, 1);
	initMaps(vorticityField, size, 0.0, 1);
	initMaps(pressureField, size, 0.0, 1);
	initMaps(divergenceField, size, 0.0, 1);
	initMaps(velField, size * 2, 0.0, 2);
	initMaps(colorField, size * 3, 0.0, 3);
	initMaps(resImage, size * 3, 0.0, 3);

	initCharMap(characteristicMap);
	initCharMap(BFECC_charMap);
	initCharMap(MM_charMap);
}

FluidSim::FluidSim(float force[2], float viscosity, float st_strength, float vort_strength, float presSteps, float projSteps, Advection advectType)
{
	this->force[0] = force[0];
	this->force[1] = force[1];
	this->visc = viscosity;
	this->st_strength = st_strength;
	this->vorticity_strength = vort_strength;

	this->advectType = advectType;

	pressureSteps = presSteps;
	projectionSteps = projSteps;

}

FluidSim::~FluidSim()
{
}

void FluidSim::initCharMap(float *&charMap)
{
	for (int j = 0; j < (nGridY + 2); j++)
	{
		for (int i = 0; i < (nGridX + 2); i++)
		{
			int ij = i + j*(nGridX + 2);

			charMap[ij * 2 + 0] = i*dXSize;
			charMap[ij * 2 + 1] = j*dYSize;

		}
	}
}

void FluidSim::sl_Advect(float *&oldField, float *&newField, size_t vecSize)
{
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int index = i + ((nGridX + 2)*j);
			float x = i*dXSize - (velField[index * 2 + 0] * dt);
			float y = j*dYSize - (velField[index * 2 + 1] * dt);

			if (x < 0.5) { x = 0.5; } if (x >(nGridX + 0.5)) { x = nGridX + 0.5; }
			int i_prev = (int)x; int i_next = i_prev + 1;
			if (y < 0.5) { y = 0.5; } if (y >(nGridY + 0.5)) { y = nGridY + 0.5; }
			int j_prev = (int)y; int j_next = j_prev + 1;

			//Bilinnear interpolation
			//Finding position inside the grid rescalled from 0-1
			float s1 = x - i_prev; float s0 = 1.0 - s1;
			float t1 = y - j_prev; float t0 = 1.0 - t1;

			int ij = index;
			int i0j0 = i_prev + (j_prev*(nGridX + 2));
			int i0j1 = i_prev + (j_next*(nGridX + 2));
			int i1j0 = i_next + (j_prev*(nGridX + 2));
			int i1j1 = i_next + (j_next*(nGridX + 2));

			
			if (vecSize == 1)
			{
				newField[ij] =	s0*(t0*oldField[i0j0] + t1*oldField[i0j1]) +
								s1*(t0*oldField[i1j0] + t1*oldField[i1j1]);
			}
			
			if (vecSize == 2)
			{
				newField[(ij * 2) + 0] = s0*(t0*oldField[i0j0 * 2 + 0] + t1*oldField[i0j1 * 2 + 0]) +
										 s1*(t0*oldField[i1j0 * 2 + 0] + t1*oldField[i1j1 * 2 + 0]);
				newField[(ij * 2) + 1] = s0*(t0*oldField[i0j0 * 2 + 1] + t1*oldField[i0j1 * 2 + 1]) +
										 s1*(t0*oldField[i1j0 * 2 + 1] + t1*oldField[i1j1 * 2 + 1]);
			}

			if (vecSize == 3)
			{
				newField[(ij * 3) + 0] = s0*(t0*oldField[i0j0 * 3 + 0] + t1*oldField[i0j1 * 3 + 0]) +
										 s1*(t0*oldField[i1j0 * 3 + 0] + t1*oldField[i1j1 * 3 + 0]);

				newField[(ij * 3) + 1] = s0*(t0*oldField[i0j0 * 3 + 1] + t1*oldField[i0j1 * 3 + 1]) +
										 s1*(t0*oldField[i1j0 * 3 + 1] + t1*oldField[i1j1 * 3 + 1]);

				newField[(ij * 3) + 2] = s0*(t0*oldField[i0j0 * 3 + 2] + t1*oldField[i0j1 * 3 + 2]) +
										 s1*(t0*oldField[i1j0 * 3 + 2] + t1*oldField[i1j1 * 3 + 2]);
			}
		}
	}
}

void FluidSim::boundaryCheck()
{
	//BOUNDARY STATES
	//Bouncing off of the Edge
	for (int j = 0; j <= nGridY+1; j++)
	{
		#pragma omp parallel for
		for (int i = 0; i <= nGridX+1; i++)
		{
			int ij = i + (j*(nGridX+2));
			if (i <= 1)
			{

				if (velField[ij * 2 + 0] < 0) { velField[ij * 2 + 0] = -velField[ij * 2 + 0]; } // LEFT MOST EDGE
			}
			else if (i >= nGridX)
			{
				if (velField[ij * 2 + 0] > 0){ velField[ij * 2 + 0] = -velField[ij * 2 + 0]; }// RIGHT MOST EDGE
			}
			if (j <= 1)
			{
				if (velField[ij * 2 + 1] < 0) { velField[ij * 2 + 1] = -velField[ij * 2 + 1]; } // TOP MOST EDGE
			}
			else if (j >= (nGridY))
			{
				if (velField[ij * 2 + 1] > 0) velField[ij * 2 + 1] = -velField[ij * 2 + 1]; // BOTTOM MOST EDGE
			}

		}
	}
}


void FluidSim::semiLangrangianAdvect(float *&oldDensity, float *&oldVelocity, float *&oldColor )
{
	
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int index = i + ((nGridX + 2)*j);
			float x = i*dXSize - (velField[index*2 + 0] * dt);
			float y = j*dYSize - (velField[index*2 + 1] * dt);

			if (x < 0.5) { x = 0.5; } if (x >(nGridX + 0.5)) { x = nGridX + 0.5; }
			int i_prev = (int)x; int i_next = i_prev + 1;
			if (y < 0.5) { y = 0.5; } if (y >(nGridY + 0.5)) { y = nGridY + 0.5; }
			int j_prev = (int)y; int j_next = j_prev + 1;

			//Bilinnear interpolation
			//Finding position inside the grid rescalled from 0-1
			float s1 = x - i_prev; float s0 = 1.0 - s1;
			float t1 = y - j_prev; float t0 = 1.0 - t1;

			int ij = index;
			int i0j0 = i_prev + (j_prev*(nGridX + 2));
			int i0j1 = i_prev + (j_next*(nGridX + 2));
			int i1j0 = i_next + (j_prev*(nGridX + 2));
			int i1j1 = i_next + (j_next*(nGridX + 2));

			//DENSITY FIELD ADVECTION
			
			densityField[ij] =	s0*(t0*oldDensity[i0j0] + t1*oldDensity[i0j1]) +
								s1*(t0*oldDensity[i1j0] + t1*oldDensity[i1j1]);
			
			//COLOR FIELD ADVECTION
			//Red 

			colorField[(ij * 3) + 0] =	s0*(t0*oldColor[(i0j0 * 3) + 0] + t1*oldColor[(i0j1 * 3) + 0]) +
									s1*(t0*oldColor[(i1j0 * 3) + 0] + t1*oldColor[(i1j1 * 3) + 0]);
		
			//Green

			colorField[(ij * 3) + 1] =	s0*(t0*oldColor[(i0j0 * 3) + 1] + t1*oldColor[(i0j1 * 3) + 1]) +
										s1*(t0*oldColor[(i1j0 * 3) + 1] + t1*oldColor[(i1j1 * 3) + 1]);

			//Blue


			colorField[(ij * 3) + 2] =	s0*(t0*oldColor[(i0j0 * 3) + 2] + t1*oldColor[(i0j1 * 3) + 2]) +
										s1*(t0*oldColor[(i1j0 * 3) + 2] + t1*oldColor[(i1j1 * 3) + 2]);
			//VELOCITY FIELD ADVECTION
			//Velocity X

			velField[(ij * 2) + 0] =	s0*(t0*oldVelocity[(i0j0 * 2) + 0] + t1*oldVelocity[(i0j1 * 2) + 0]) +
										s1*(t0*oldVelocity[(i1j0 * 2) + 0] + t1*oldVelocity[(i1j1 * 2) + 0]);
			//Velocity Y

			velField[(ij * 2) + 1] =	s0*(t0*oldVelocity[(i0j0 * 2) + 1] + t1*oldVelocity[(i0j1 * 2) + 1]) +
										s1*(t0*oldVelocity[(i1j0 * 2) + 1] + t1*oldVelocity[(i1j1 * 2) + 1]);

			//BOUNDARY STATES
			//Bouncing off of the Edge
			if (i == 1)
			{

				if (velField[ij * 2 + 0] < 0) velField[ij * 2 + 0] = -velField[ij * 2 + 0]; // LEFT MOST EDGE
			}
			else if (i == nGridX)
			{
				if (velField[ij * 2 + 0] > 0) velField[ij * 2 + 0] = -velField[ij * 2 + 0]; // RIGHT MOST EDGE
			}
			if (j == 1)
			{
				if (velField[ij * 2 + 1] < 0) velField[ij * 2 + 1] = -velField[ij * 2 + 1]; // TOP MOST EDGE
			}
			else if (j == nGridY)
			{
				if (velField[ij * 2 + 1] > 0) velField[ij * 2 + 1] = -velField[ij * 2 + 1]; // BOTTOM MOST EDGE
			}

		}
	}
}

void FluidSim::interpVelocity(float x, float y, float &resVelx, float &resVely, float *&velocityField)
{
//	if (x < -(nGridX + 0.5)) { x = 0.5; } if (x > 2 * (nGridX + 0.5)) { x = (nGridX + 0.5);}
//	if (x < 0.5) { x = 0.5 - x; } if (x >(nGridX + 0.5)) { x = 2 * (nGridX + 0.5) - x; }	
	if (x < 0.5) { x = 0.5; } if (x >(nGridX + 0.5)) { x = (nGridX + 0.5) ; }
	int i_prev = (int)x; int i_next = i_prev + 1;
//	if (y < -(nGridY + 0.5)) { y = 0.5; } if (y > 2 * (nGridY + 0.5)) { y = (nGridY + 0.5); }
//	if (y < 0.5) { y = 0.5 - y; } if (y >(nGridY + 0.5)) { y = 2*(nGridY + 0.5)-y; }
	if (y < 0.5) { y = 0.5; } if (y >(nGridY + 0.5)) { y = (nGridY + 0.5); }
	int j_prev = (int)y; int j_next = j_prev + 1;

	//Bilinnear interpolation
	//Finding position inside the grid rescalled from 0-1
	float s1 = x - i_prev; float s0 = 1.0 - s1;
	float t1 = y - j_prev; float t0 = 1.0 - t1;

	
	int i0j0 = i_prev + (j_prev*(nGridX + 2));
	int i0j1 = i_prev + (j_next*(nGridX + 2));
	int i1j0 = i_next + (j_prev*(nGridX + 2));
	int i1j1 = i_next + (j_next*(nGridX + 2));

	//Velocity X
	resVelx = s0*(t0*velocityField[(i0j0 * 2) + 0] + t1*velocityField[(i0j1 * 2) + 0]) +
			  s1*(t0*velocityField[(i1j0 * 2) + 0] + t1*velocityField[(i1j1 * 2) + 0]);
	//Velocity Y
	resVely = s0*(t0*velocityField[(i0j0 * 2) + 1] + t1*velocityField[(i0j1 * 2) + 1]) +
			  s1*(t0*velocityField[(i1j0 * 2) + 1] + t1*velocityField[(i1j1 * 2) + 1]);


}

void FluidSim::MM_AdvectCharMap()
{
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int ij = i + (j*(nGridX + 2));
			float x = characteristicMap[ij * 2 + 0];
			float y = characteristicMap[ij * 2 + 1];

			float x_forward = x - velField[ij * 2 + 0] * dt;
			float y_forward = y - velField[ij * 2 + 1] * dt;
			float velx_forward, vely_forward;
			velx_forward = vely_forward = 0;
			interpVelocity(x_forward, y_forward, velx_forward, vely_forward, velField);
			float x_backward = x_forward + velx_forward * dt;
			float y_backward = y_forward + vely_forward * dt;

			float error_x = 0.5f * (x - x_backward);
			float error_y = 0.5f * (y - y_backward);

			float x_mm = x_forward + error_x;
			float y_mm = y_forward + error_y;
//			if (x_mm < -(nGridX + 0.5)) { x_mm = 0.5; } if (x_mm > 2 * (nGridX + 0.5)) { x_mm = (nGridX + 0.5); }
//			if (x_mm < 0.5) { x_mm = 0.5 - x; } if (x_mm >(nGridX + 0.5)) { x_mm = 2*(nGridX + 0.5) - x ; }
			
//			if (y_mm < -(nGridY + 0.5)) { y_mm = 0.5; } if (y_mm > 2 * (nGridY + 0.5)) { y_mm = (nGridY + 0.5); }
//			if (y_mm < 0.5) { y_mm = 0.5 - y; } if (y_mm >(nGridY + 0.5)) { y_mm = 2*(nGridY + 0.5) - y ; }
			if (x_mm < 0.5) { x_mm = 0.5; } if (x_mm >(nGridX + 0.5)) { x_mm = nGridX + 0.5; }
			if (y_mm < 0.5) { y_mm = 0.5; } if (y_mm >(nGridY + 0.5)) { y_mm = (nGridY + 0.5); }

			MM_charMap[ij * 2 + 0] = x_mm;
			MM_charMap[ij * 2 + 1] = y_mm;
/*			if (ij == 513)
			{
				std::cout << "i: " << i << " j: " << j << std::endl;
				std::cout << "x: " << x << " y: " << y << std::endl;
				std::cout << "velX: " << velField[ij * 2 + 0] << " velY: " << velField[ij * 2 + 1] << std::endl;
				std::cout << "x_foreward: " << x_forward << " y_forward: " << y_forward << std::endl;
				std::cout << "x_backward: " << x_backward << " y_backward: " << y_backward << std::endl;
				std::cout << "error_x: " << error_x << " error_y: " << error_y << std::endl;
				std::cout << "x_mm: " << x_mm << " y_mm: " << y_mm << std::endl;
				
			}
*/
		}
	}

}

void FluidSim::BFECC_AdvectCharMap()
{
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int ij = i + j*(nGridX + 2);
			float x = characteristicMap[ij * 2 + 0];
			float y = characteristicMap[ij * 2 + 1];
			float x_forward = x - velField[ij * 2 + 0] * dt;
			float y_forward = y - velField[ij * 2 + 1] * dt;
			float velx_forward, vely_forward;
			velx_forward = vely_forward = 0;
			interpVelocity(x_forward, y_forward, velx_forward, vely_forward, velField);
			float x_backward = x_forward + velx_forward * dt;
			float y_backward = y_forward + vely_forward * dt;
			
			float error_x = 0.5f * (x - x_backward);
			float error_y = 0.5f * (y - y_backward);

			float x_bfe = x + error_x;
			float y_bfe = y + error_y;

			float velx_bfe, vely_bfe;
			velx_bfe = vely_bfe = 0;
			interpVelocity(x_bfe, y_bfe, velx_bfe, vely_bfe, velField);
			float x_bfecc = x_bfe - velx_bfe * dt;
			float y_bfecc = y_bfe - vely_bfe * dt;
			
			BFECC_charMap[ij * 2 + 0] = x_bfecc;
			BFECC_charMap[ij * 2 + 1] = y_bfecc;
		}
	}

}

void FluidSim::charMap_Advect(float *&charMap, float *&oldField, float *&newField, size_t vecSize)
{
	for (int j = 1; j <= nGridY; j++)
	{

		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int ij = i + j*(nGridX + 2);

			float x = charMap[ij * 2 + 0];
			float y = charMap[ij * 2 + 1];
			//std::cout << "x : " << x << " y: " << y << std::endl;
//			if (x < -(nGridX + 0.5)) { x = 0.5; } if (x > 2 * (nGridX + 0.5)) { x = (nGridX + 0.5); }
//			if (x < 0.5) { x = 0.5 - x; } if (x >(nGridX + 0.5)) { x = 2*(nGridX + 0.5) - x ; }
			if (x < 0.5) { x = 0.5; } if (x >(nGridX + 0.5)) { x = (nGridX + 0.5); }
			int i_prev = (int)x; int i_next = i_prev + 1;
//			if (y < -(nGridY + 0.5)) { y = 0.5; } if (y > 2 * (nGridY + 0.5)) { y = (nGridY + 0.5); }
//			if (y < 0.5) { y = 0.5 - y; } if (y >(nGridY + 0.5)) { y = 2*(nGridY + 0.5) - x ; }
			if (y < 0.5) { y = 0.5; } if (y >(nGridY + 0.5)) { y = (nGridY + 0.5); }
			int j_prev = (int)y; int j_next = j_prev + 1;

			//Bilinnear interpolation
			//Finding position inside the grid rescalled from 0-1
			float s1 = x - i_prev; float s0 = 1.0 - s1;
			float t1 = y - j_prev; float t0 = 1.0 - t1;

			
			int i0j0 = i_prev + (j_prev*(nGridX + 2));
			int i0j1 = i_prev + (j_next*(nGridX + 2));
			int i1j0 = i_next + (j_prev*(nGridX + 2));
			int i1j1 = i_next + (j_next*(nGridX + 2));


			if (vecSize == 1)
			{
				newField[ij] = s0*(t0*oldField[i0j0] + t1*oldField[i0j1]) +
							   s1*(t0*oldField[i1j0] + t1*oldField[i1j1]);
			}

			if (vecSize == 2)
			{
				newField[(ij * 2) + 0] = s0*(t0*oldField[i0j0 * 2 + 0] + t1*oldField[i0j1 * 2 + 0]) +
					s1*(t0*oldField[i1j0 * 2 + 0] + t1*oldField[i1j1 * 2 + 0]);
				newField[(ij * 2) + 1] = s0*(t0*oldField[i0j0 * 2 + 1] + t1*oldField[i0j1 * 2 + 1]) +
					s1*(t0*oldField[i1j0 * 2 + 1] + t1*oldField[i1j1 * 2 + 1]);
			}

			if (vecSize == 3)
			{
				newField[(ij * 3) + 0] = s0*(t0*oldField[i0j0 * 3 + 0] + t1*oldField[i0j1 * 3 + 0]) +
					s1*(t0*oldField[i1j0 * 3 + 0] + t1*oldField[i1j1 * 3 + 0]);

				newField[(ij * 3) + 1] = s0*(t0*oldField[i0j0 * 3 + 1] + t1*oldField[i0j1 * 3 + 1]) +
					s1*(t0*oldField[i1j0 * 3 + 1] + t1*oldField[i1j1 * 3 + 1]);

				newField[(ij * 3) + 2] = s0*(t0*oldField[i0j0 * 3 + 2] + t1*oldField[i0j1 * 3 + 2]) +
					s1*(t0*oldField[i1j0 * 3 + 2] + t1*oldField[i1j1 * 3 + 2]);
			}
		}
	}

}





void FluidSim::advectionStep()
{
	int size = (nGridX + 2)*(nGridY + 2);
	float* tempDensity = new float[size];
	float* tempVelocity = new float[size*2];
	float* tempColor = new float[size*3];

	memcpy(tempDensity, densityField, size * sizeof(float));
	memcpy(tempVelocity, velField, size * 2 * sizeof(float));
	memcpy(tempColor, colorField, size * 3 * sizeof(float));

	
	if (advectType == Advection::SL)
	{
		sl_Advect(tempDensity, densityField, 1);
		sl_Advect(tempColor, colorField, 3);
		sl_Advect(tempVelocity, velField, 2);
		
	}

	else if (advectType == Advection::BFECC)
	{
		BFECC_AdvectCharMap();
		charMap_Advect(BFECC_charMap, tempDensity, densityField, 1);
		charMap_Advect(BFECC_charMap, tempColor, colorField, 3);
		charMap_Advect(BFECC_charMap, tempVelocity, velField, 2);
	}
	
	else if (advectType == Advection::MM)
	{
		MM_AdvectCharMap();
		charMap_Advect(MM_charMap, tempDensity, densityField, 1);
		charMap_Advect(MM_charMap, tempColor, colorField, 3);
		charMap_Advect(MM_charMap, tempVelocity, velField, 2);

	}
//	boundaryCheck();
	//semiLangrangianAdvect(tempDensity, tempVelocity, tempColor);

	delete[] tempDensity;
	delete[] tempVelocity; 
	delete[] tempColor;

}


void FluidSim::viscosityCalc(int index_i, int index_j, float &xVel, float &yVel)
{
	int ij = index_i + (index_j*(nGridX + 2));
	int i1j = (index_i + 1) + (index_j*(nGridX + 2));
	int i0j = (index_i - 1) + (index_j*(nGridX + 2));
	int ij1 = index_i + ((index_j + 1)*(nGridX + 2));
	int ij0 = index_i + ((index_j - 1)*(nGridX + 2));

	float e = exp(-4.0f * visc*dt);
	float xField = velField[i0j * 2 + 0] + velField[i1j * 2 + 0] + velField[ij0 * 2 + 0] + velField[ij1 * 2 + 0];
	float yField = velField[i0j * 2 + 1] + velField[i1j * 2 + 1] + velField[ij0 * 2 + 1] + velField[ij1 * 2 + 1];
	xVel = e*(velField[ij * 2 + 0] + (visc*dt*xField) / (dXSize*dXSize));
	yVel = e*(velField[ij * 2 + 1] + (visc*dt*yField) / (dXSize*dXSize));
	//if (abs(xField) > 10) std::cout << xField << std::endl;
}


void FluidSim::SurfaceTensionCalc(int index_i, int index_j, float &xVel, float &yVel)
{
	int ij = index_i + (index_j*(nGridX + 2));
	int i1j = (index_i + 1) + (index_j*(nGridX + 2));
	int i0j = (index_i - 1) + (index_j*(nGridX + 2));
	int ij1 = index_i + ((index_j + 1)*(nGridX + 2));
	int ij0 = index_i + ((index_j - 1)*(nGridX + 2));
	
	int i1j1 = (index_i + 1) + (index_j + 1)*(nGridX + 2);
	int i0j0 = (index_i - 1) + (index_j - 1)*(nGridX + 2);
	int i1j0 = (index_i + 1) + (index_j - 1)*(nGridX + 2);
	int i0j1 = (index_i - 1) + (index_j + 1)*(nGridX + 2);


	float grad_X = (densityField[i1j] - densityField[i0j]) / (2.0*dXSize);
	float grad_Y = (densityField[ij1] - densityField[ij0]) / (2.0*dYSize);
	float normalizer = glm::sqrt(grad_X*grad_X + grad_Y*grad_Y);

	if (normalizer >= 1e-10)
	{
		float mxx = (densityField[i1j] + densityField[i0j] - 2 * densityField[ij]) / (dXSize * dXSize);
		float myy = (densityField[i1j] + densityField[i0j] - 2 * densityField[ij]) / (dYSize * dYSize);
		float mxy = (densityField[i1j1] + densityField[i0j0] - densityField[i1j0] - densityField[i0j1]) / (dXSize * dYSize * 4);


		float Kij = (mxx + myy - (grad_X*grad_X*mxx) - (grad_Y*grad_Y*myy) + (2 * grad_X*grad_Y*mxy)) / normalizer;

		grad_X = grad_X / normalizer; grad_Y = grad_Y / normalizer;
		xVel = Kij*grad_X*st_strength*dt;
		yVel = Kij*grad_Y*st_strength*dt;
		//if (fabs(xVel) > 0) std::cout << xVel << std::endl;
	}
}

void FluidSim::vorticityConfinementCalc(int index_i, int index_j, float &xVel, float &yVel)
{
	int ij = index_i + index_j*(nGridX + 2);
	int i1j = (index_i + 1) + index_j*(nGridX + 2);
	int i0j = (index_i - 1) + index_j*(nGridX + 2);
	int ij0 = index_i + (index_j - 1)*(nGridX + 2);
	int ij1 = index_i + (index_j + 1)*(nGridX + 2);

	float gradX = (abs(vorticityField[i1j]) - abs(vorticityField[i0j])) / (2.0f * dXSize);
	float gradY = (abs(vorticityField[ij1]) - abs(vorticityField[ij0])) / (2.0f * dXSize);

	float normalizer = glm::sqrt(gradX*gradX + gradY*gradY);
	if (normalizer > 1e10)
	{
		gradX = gradX / normalizer;
		gradY = gradY / normalizer;

		xVel = gradX*vorticityField[ij] * vorticity_strength*dXSize * dt;
		yVel = gradY*vorticityField[ij] * vorticity_strength*dYSize * dt;
	}
	
}

void FluidSim::addingSources(float *&sourceDensity)
{
	int size = (nGridX + 2)*(nGridY + 2);
	#pragma omp parallel for
	for (int i = 0; i < size; i++)
	{
		float temp = sourceDensity[i]*dt;
		float temp2 = densityField[i];
		densityField[i] = temp2 + temp;
	}
}

void FluidSim::addingForces()
{
	//Bouyancy
	int size = (nGridX + 2)*(nGridY + 2);
	
	for (int j = 1; j <= (nGridY) ; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= (nGridX); i++)
		{
			int index = i + j * (nGridX + 2);
			
			if (vorticity_strength != 0)
			{
				float xVel, yVel;
				xVel = yVel = 0;
				vorticityConfinementCalc(i, j, xVel, yVel);
				velField[(index * 2) + 0] += xVel;
				velField[(index * 2) + 1] += yVel;
			}
			if (st_strength != 0)
			{
				float xVel, yVel;
				xVel = yVel = 0;
				SurfaceTensionCalc(i, j, xVel, yVel);
				velField[(index * 2) + 0] += xVel;
				velField[(index * 2) + 1] += yVel;
			}

			if (visc != 0)
			{
				float xVel, yVel;
				xVel = yVel = 0;
				viscosityCalc(i, j, xVel, yVel);
				velField[(index * 2) + 0] += xVel;
				velField[(index * 2) + 1] += yVel;
			}
			
			float xF = densityField[index] * force[0];
			float yF = densityField[index] * force[1];

			velField[(index * 2) + 1] += yF*dt;
			velField[(index * 2) + 0] += xF*dt;

			
		}
	}
}

void FluidSim::calcVorticity()
{
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int ij = i + j*(nGridX + 2);
			int i0j = (i - 1) + j*(nGridX + 2);
			int i1j = (i + 1) + j*(nGridX + 2);
			int ij0 = i + (j - 1)*(nGridX + 2);
			int ij1 = i + (j + 1)*(nGridX + 2);

			vorticityField[ij] = ((velField[ij1 * 2 + 0] - velField[ij0 * 2 + 0]) / (2.0f*dYSize)) - ((velField[i1j*2 + 1] - velField[i0j*2 + 1])/(2.0f*dXSize));

		}
	}
}

void FluidSim::calcDivergence()
{
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int ij = i + (j*(nGridX + 2));
			int i1j = (i + 1) + (j*(nGridX + 2));
			int i0j = (i - 1) + (j*(nGridX + 2));

			int ij1 = i + ((j + 1)*(nGridX + 2));
			int ij0 = i + ((j - 1)*(nGridX + 2));

			float divX = 0.5 * (velField[i1j * 2] - velField[i0j * 2])*(1.0f / dXSize); //X velocity Field

			float divY = 0.5 * (velField[(ij1 * 2) + 1] - velField[(ij0 * 2) + 1])*(1.0f / dYSize);//Y velocity Field
			divergenceField[ij] = divX + divY;
		}
	}
}

void FluidSim::projection()
{
	calcDivergence();
	initMaps(pressureField,(nGridY+2)*(nGridX+2), 0.0, 1);
	//Pressure Step
	for (int iter = 0; iter < pressureSteps; iter++)
	{
		for (int j = 1; j <= nGridY; j++)
		{
			for (int i = 1; i <= nGridX; i++)
			{
				int ij = i + (j*(nGridX + 2));
				int i1j = (i + 1) + (j*(nGridX + 2));
				int i0j = (i - 1) + (j*(nGridX + 2));

				int ij1 = i + ((j + 1)*(nGridX + 2));
				int ij0 = i + ((j - 1)*(nGridX + 2));

				float divergence = divergenceField[ij];

				float avgPressure = 0.25f * (pressureField[i1j] + pressureField[i0j] + pressureField[ij1] + pressureField[ij0]);
				pressureField[ij] = avgPressure - (0.25f * divergence * dXSize * dXSize);
			}
		}
	}

	//Velocity Step
	for (int j = 1; j <= nGridY; j++)
	{
		#pragma omp parallel for
		for (int i = 1; i <= nGridX; i++)
		{
			int ij = i + (j*(nGridX + 2));
			int i1j = (i + 1) + (j*(nGridX + 2));
			int i0j = (i - 1) + (j*(nGridX + 2));

			int ij1 = i + ((j + 1)*(nGridX + 2));
			int ij0 = i + ((j - 1)*(nGridX + 2));

			float pressure_gradiant_X = 0.5f * (1.0f / dXSize) * (pressureField[i1j] - pressureField[i0j]);
			float pressure_gradiant_Y = 0.5f * (1.0f / dYSize) * (pressureField[ij1] - pressureField[ij0]);

			velField[(ij * 2) + 0] = velField[(ij * 2) + 0] - pressure_gradiant_X;
			velField[(ij * 2) + 1] = velField[(ij * 2) + 1] - pressure_gradiant_Y;
		}
	}
}

void FluidSim::boundaryConditions(float *&sourceObstruction)
{
	for (int j = 1; j <= nGridY + 1; j++)
	{
#pragma omp parallel for
		for (int i = 1; i <= nGridX + 1; i++)
		{
			int index = i + (j*(nGridX + 2));

						if (i <= 1) {
							if (velField[index * 2 + 0] < 0.0f) velField[index * 2 + 0] = -velField[index * 2 + 0];
						}
						if (i >= nGridX)
						{
							if (velField[index * 2 + 0] > 0.0f) velField[index * 2 + 0] = -velField[index * 2 + 0];
						}
						if (j <= 1)
						{
							if (velField[index * 2 + 1] < 0.0f) velField[index * 2 + 1] = -velField[index * 2 + 1];
						}
						if (j >= nGridY)
						{
							if (velField[index * 2 + 1] > 0.0f) velField[index * 2 + 1] = -velField[index * 2 + 1];
						}
			
//			if (i <= 1 || i >= nGridX) velField[index * 2 + 0] = 0;
//			if (j <= 1 || j >= nGridY) velField[index * 2 + 1] = 0;
		
			velField[index * 2 + 0] *= sourceObstruction[index];
			velField[index * 2 + 1] *= sourceObstruction[index];
		}
	}
}


void FluidSim::simulate(float *&inputImage, float *&sourceDensity, float *&sourceObstruction, float dt)
{
	this->colorField = inputImage;
	this->dt = dt;
	
	advectionStep();
	calcVorticity();
	addingSources(sourceDensity);
	addingForces();
	
	for (int i = 0; i < projectionSteps; i++)
	{
		boundaryConditions(sourceObstruction);
		projection();
	}

	resImage = colorField;
}