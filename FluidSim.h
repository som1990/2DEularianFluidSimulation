/*
Author: Soumitra Goswami
Date: 4/8/2018

Description: This simulation class simulates a 2D Fluid Simulation giving you the option of 3 advection Schemes:
                 - Semi Lagrangian Advection
	             - Backward Forward Error Compensation and Correction Method(BFECC)
	             - Modified MacCormack Advection
             The program also generates a velocity field, pressure field, Divergence Field, Density field and color field
			 that can be visualized by calling the public fields.
			 The program supports usage of Characteristic Maps, Incompressability, Vorticity Confinement and Obstruction Maps.
			 

*/

#ifndef FLUIDSIM_H
#define FLUIDSIM_H
#include <iostream>
#include <glm/glm.hpp>

class FluidSim
{
public:
// Accessable Variables
	float *resImage;
	float *densityField, *velField, *colorField, *pressureField, *divergenceField; // fields generated
	enum Advection { BFECC, SL, MM };//Supports 3 advection schemes
	Advection advectType;

// Constructors	
	FluidSim(float force[2], float viscocity, float ST_strength, float vort_strength, float presSteps, float projSteps, Advection advectionType );
	~FluidSim();
	
// Accessable Functions	
    /* void initFields(int,int,float,float)
    //Description: Initializes the all the grid fields to the sizes provided
    //Input: 
    //      @int nGridX : Num grids in X direction
    //      @int nGridY : Num grids in Y direction
	//      @float dXSize : Horizontal grid cell size
    //      @float dYSize : Vertical grid cell size
    //Returns: null (Just sets the internal fields in the class)
	*/
	void initFields(int nGridsX, int nGridsY, float dXSize, float dYSize);
	
	
    /* void simulate(float* , float* , float*, float)
    //Description: This is where our main simulation step takes place.
	//Input:
	//      @float* inputImage : The color field to advect our fluid simulation. 
	//		@float* sourceDensity: The density field used to advect our simulation.
	//      @float* sourceObstruction: The obstruction field to use to act as our boundry.
    //      @float dt : Time step to use to simulate our fluid.
	//Returns: null (Updates the fields based on the timestep)
    */
	void simulate(float *&inputImage, float *&sourceDensity, float *&sourceObstruction, float dt);
	
	

private:
	//Inaccessable Variables
	float force[2], dXSize, dYSize, visc, st_strength, vorticity_strength, dt;
	float *characteristicMap, *BFECC_charMap, *MM_charMap, *vorticityField;
	int pressureSteps, projectionSteps;
	int nGridX, nGridY;

	//Inaccessable functions
	
	/* void AdvectionStep()
	// Description: Sends the data to different advection schemes based on user input 
	*/
	void advectionStep();
	
	/*void addingSources(float* )
	// Description: Takes the user painted source map and adds it to the current density field.
	// Input: @float* sourceDensity - the painted source field to add to the current density.
	*/
	void addingSources(float *&sourceDensity);
	
	/*void AddingForces()
	//Description: Adds the calculated external forces (Viscosity , Vorticity Confinement , bouyancy) and updates
	// the velocity field
	*/
	void addingForces();
	
	/*void calcDivergence()
	// Description: Calculates the divergence of the density field which will later be used to make the fluid 
	//				Incompressable in the Pressure Projection Step. This compensation creates the curly nature of the fluid.
	//				
	*/
	void calcDivergence();
	
	/* void calcVorticity()
	// Description: Calculates the vorticity ( the curl of the simulation) and then can be made more swirly or less swirly based
	//              on user input.
	*/
	void calcVorticity();
	
	/* void projection()
	// Description: Updates the velocity by calculating the Gauss Seidel pressure projection and compensating it by subtracting divergence.
	//              This method is really slow. So better but more complicated methods could be utilized to make this step much faster 
	                eg:(Conjugate Gradient Method or Multigrid).
	*/
	void projection();
	
	/* void viscocity(int, int, float, float) 
	// Description: Updates the velocity by calculating the viscocity of the model. Thicker fluids move slower, 
	//              Thinner fluids move faster
	// Input: @ int index_i : The current horizontal grid index
	//        @ int index_j : The current vertical grid index
	//        @ float xVel  : Horizontal velocity to update 
	//        @ float yVel  : Vertical velocity to update
	//
	*/
	void viscosityCalc(int index_i, int index_j, float &xVel, float &yVel);
	
	/* void surfaceTensionCalc(int, int, float, float) 
	// Description: Updates the velocity by calculating the surfaceTension of the fluid. Higher surfacetension the 
	                more the fluid will stick together as a blob. Less surfaceTension makes the fluid more wispy and break apart easily
	//              
	// Input: @ int index_i : The current horizontal grid index
	//        @ int index_j : The current vertical grid index
	//        @ float xVel  : Horizontal velocity to update 
	//        @ float yVel  : Vertical velocity to update
	//
	*/
	void SurfaceTensionCalc(int index_i, int index_j, float &xVel, float &yVel);
	
	/* void vorticityConfinementCalc(int, int, float, float) 
	// Description: Calculates vorticity and allows the user to increase/decrease the swirliness of the fluid 
	// Input: @ int index_i : The current horizontal grid index
	//        @ int index_j : The current vertical grid index
	//        @ float xVel  : Horizontal velocity to update 
	//        @ float yVel  : Vertical velocity to update
	//
	*/
	void vorticityConfinementCalc(int index_i, int index_j , float &xVel, float &yVel);
	
	/* void initCharMap(float* )
	// Description: Initializes the characteristic Map by filling the grid with grid position (i,j).
	// Input : @float* charMap - The charasteristic Map to populate.
	*/
	void initCharMap(float *&charMap);
	
	/* void BFECC_AdvectCharMap()
	// Description: It utilizes the characteristic map for advection. Using the Backward Forward Error Compensation & Correction Method.
	//              It occers in 5 steps:
					1) We do a Semi Lagrangian Advection forward 
						f_forward(x) =  f0(x - u(x)dt)
					2) We do a Semi Lagrangian Backward
					    f_backward(x) = f_forward(x + u(x)dt)
					3) We calculate the error
                       err(x) = 0.5 * (f0(x) - f_backward(x))
                    4) We compensate for the error
                       f_bfe(x) = f0(x) + err(x)
                    5) We then do an semi lagrangian advection on the compensated functions to get our final advected result
                       f_bfecc(x) = f_bfe(x-u(x)dt)
					    
	// 
	*/
	void BFECC_AdvectCharMap();
	
	/* void MM_AdvectCharMap()
	// Description: It utilizes the characteristic map for advection. Using the Modified MacCormack advection scheme.
	//              Similar to BFECC but change in step 4 . It occers in 4 steps:
					1) We do a Semi Lagrangian Advection forward 
						f_forward(x) =  f0(x - u(x)dt)
					2) We do a Semi Lagrangian Backward
					    f_backward(x) = f_forward(x + u(x)dt)
					3) We calculate the error
                       err(x) = 0.5 * (f0(x) - f_backward(x))
                    4) We compensate for the error adding it to the forward advection
                       f_mm(x) = f_forward(x) + err(x)
					    
	// 
	*/
	void MM_AdvectCharMap();
	
	/* interpVelocity(float,float,float,float,float)
	// Description: Does a bilinear interpolation to update the velocity at a given grid point (x,y)
	// Input: @float x : horizontal spacial position within the grid
              @float y : vertical spacial position within the grid
			  @float resVelx : updated calculated horizontal velocity at the position(x,y)
			  @float resVely : updated calculated vertical velocity at the position(x,y)
			  @float* velocityField: Gridded Velocity field used to sample the values 
	*/
	void interpVelocity(float x, float y, float &resVelx, float &resVely, float *&velocityField);
	
	/*charMap_Advect( float*, float*, float*, int)
	// Description: Using the advected characteristic Maps ( these store advected positional values) to bilinear interpolate 
	                the fields based on their dimension. 
	// Input: @float* charMap - The already advected charmap based on the methods (In our case on BFECC and MM ).
	//        @float* oldField - The field used to sample values
	          @float* newField - The field to store the interpolated values 
              @size_t vecSize - Dimension of the field ( 1D - for density, 2D - for velocity)			  
	*/
	void charMap_Advect(float *&charMap, float *&oldField, float *&newField, size_t vecSize);
	
	/*void SemiLagrangianAdvection( float*, float* , float*)
	// Description: Advects based on back-tracing the velocity values. Read the stable fluids paper for more details
	                http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf
					f(x) = x - u(x)dt;
					In this function I didn't separate each step and advected all the fields at once.
	   Input: @float* oldDensity - the density field to advect
	          @float* oldVelocity - the velocity field to advect
			  @float* oldColor - the color field to advect
	*/
	
	void semiLangrangianAdvect(float *&oldDensity, float *&oldVelocity, float *&oldColor);
	
	/*void sl_Advect( float*, float* , size_t)
	// Description: Advects based on back-tracing the velocity values. Read the stable fluids paper for more details
	                http://www.dgp.toronto.edu/people/stam/reality/Research/pdf/ns.pdf
					f(x) = x - u(x)dt;
					In this function I generalized the sl advection step based on field dimensionality. Makes it more 
					flexible to use.
	   Input: @float* oldField - the field used to sample values
	          @float* newField - the field used to store advected interpolated values.
			  @size_t vecSize - the dimension of the field.
			  
			  
	*/
	void sl_Advect(float *&oldField, float *&newField, size_t vecSize);
	
	/*void boundaryCheck()
		Description: Checks the edge boundaries of the image to bounce the fluid off the edge
	*/
	void boundaryCheck();
	
	
	/*void boundaryConditions(float*)
	Description: Takes the obstruction map and sets the velocity at those locations to zero. Just like boundaryCheck it 
	             checks for edge cases and changes the velocity direction to bounce the fluid off the edge.
	*/
	void boundaryConditions(float *&sourceObstruction);
};

#endif

