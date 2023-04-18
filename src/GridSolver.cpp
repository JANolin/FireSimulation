#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <utility>
#include <stdlib.h>

#include "GridSolver.h"

// ###############################
// TO IMPLEMENT
// ###############################
// DONE - Black Body Radiation
// DONE - Add Scalar Field for new 
// DONE - Gravity Forces
// DONE - Expansion forces
// DONE - Combustion
// DONE - Add sliders for remaining parameters
// DONE - Add fuel implementation
// 
// 
// DONE - Duplicate all diffusion parameters, etc for all density field
// DONE - Add the sliders
// DONE - Create lookup table for color
// 
// 
// DONE -MAYBE ADD VORTICITY CONFINEMENT
// 
// 
// 
// DONE - Viz Exhaust Smoke
// DONE - Test everything
// Add slider for vorcity
// Tune everything
// Tune Black Body Radiation
// ###############################




using std::vector;
using std::swap;
using std::deque;
using std::shared_ptr;
using std::make_shared;


const double pi = 3.1415926535897932385;
const double light_spd = 2.99792458e8;
const double planck_cst = 6.62606957e-34;
const double kBoltzmann = 1.380658e-23;
const double wavelenghtRed = 630e-7;
const double wavelenghtGreen = 532e-7;
const double wavelenghtBlue = 465e-7;
const int blueCoreCenter = cols / 2 + 1;

float nv = 1.0f, kappa = 0.05f; // diffuse rate. nv is for velocity field and kapp is for heat transfer.
float kappa_fuel = 0.1f; 
float kappa_exhaust = 0.1f;
float dissipation = 0.05f;
float dissipation_fuel = 0.01f;
float dissipation_exhaust = 0.05f;

float dt = 0.01f; // simulation step size
float simTime;
float dx = 1; // size of a cell in simulation
float velScale = 10.0f; // scaling for display
float buoyancy = 8.f; // buoyancy force coefficient
float refinementThreshold = 20.0f;
int iterations = 20; // number of iterations for iterative solvers (Gauss–Seidel, Conjugate Gradients, etc.)
float gravity = 0.5f;
float vorticityStrength = 1.f;
float reactionTemperature = 30;
float expansion = 0.05;
float pctFuelBurn = 0.05;
float stochiometricMixture = 1;
float combustionHeat = 10;
float mapIntensity = 100;
int vizType = 0;

bool snake_flag = false;
/*
 * Quantities in the grid.
 * Each has 2 entries (e.g., vx[0] and vy[1] for storing the current/next state)
 * But you had better use pointers like cur_vx, next_vx, etc. defined below.
 * Note that rows are horizontal and columns are vertical.
 * e.g., cur_vx[3][4] is the x component of the velocity at position (4, 3).
 * Also note that the first rows/columns and last rows/columns are for the boundary.
 */
// velocity
double vx[2][rows + 2][cols + 2];
double vy[2][rows + 2][cols + 2];

// temperature
double tp[2][rows + 2][cols + 2];

// fuel gases
double fg[2][rows + 2][cols + 2];
// exhaust gases
double eg[2][rows + 2][cols + 2];

// pointers for the corresponding quantities in the current and the next states
double (*cur_vx)[cols + 2] = vx[0], (*next_vx)[cols + 2] = vx[1];
double (*cur_vy)[cols + 2] = vy[0], (*next_vy)[cols + 2] = vy[1];

double(*cur_tp)[cols + 2] = tp[0], (*next_tp)[cols + 2] = tp[1];
double(*cur_fg)[cols + 2] = fg[0], (*next_fg)[cols + 2] = fg[1];
double(*cur_eg)[cols + 2] = eg[0], (*next_eg)[cols + 2] = eg[1];

// flow source
double vxSrc[rows + 2][cols + 2];
double vySrc[rows + 2][cols + 2];

// heat source
double tpSrc[rows + 2][cols + 2];

// fuel source
double fgSrc[rows + 2][cols + 2];

// exhaust source
double egSrc[rows + 2][cols + 2];

//filament
deque<shared_ptr<vector<Vertex>>> filaments;
deque<double> ages;
float maxAge = 10;

/*
 * Flow type:
 *	horizontal: horizontal flow
 *	vertical: vertical flow
 *	other: quantities like temperature, density, etc.
 */
enum FlowType { horizontal, vertical, other };
void setBnd(double a[rows + 2][cols + 2], FlowType flowType);

Vertex gridVertices[rows + 1][cols + 1];
GLuint gridIndices[6 * rows * cols];
Vertex velVertices[rows][cols][2];

void initGrid() {
	memset(vx, 0, sizeof(vx));
	memset(vy, 0, sizeof(vx));
	memset(tp, 0, sizeof(vx));
	memset(fg, 0, sizeof(vx));
	memset(eg, 0, sizeof(vx));
	memset(vxSrc, 0, sizeof(vxSrc));
	memset(vySrc, 0, sizeof(vySrc));
	memset(tpSrc, 0, sizeof(tpSrc));
	memset(fgSrc, 0, sizeof(tpSrc));
	memset(egSrc, 0, sizeof(tpSrc));
	filaments.clear();
	ages.clear();
	cur_vx = vx[0];
	cur_vy = vy[0];
	cur_tp = tp[0];
	cur_fg = fg[0];
	cur_eg = eg[0];
	next_vx = vx[1];
	next_vy = vy[1];
	next_tp = tp[1];
	next_fg = fg[1];
	next_eg = eg[1];


	snake_flag = false;

	GLuint indices[rows + 1][cols + 1];
	GLuint idx = 0;
	for (size_t i = 0; i <= rows; ++i) {
		for (size_t j = 0; j <= cols; ++j) {
			gridVertices[i][j] = Vertex((float)j * cellSize, (float)i * cellSize, 0.f, 0.f, 0.f, 1.f);
			indices[i][j] = idx++;
		}
	}

	size_t k = 0;
	for (size_t i = 0; i < rows; ++i) {
		for (size_t j = 0; j < cols; ++j) {
			gridIndices[k++] = indices[i][j];
			gridIndices[k++] = indices[i][j + 1];
			gridIndices[k++] = indices[i + 1][j + 1];
			gridIndices[k++] = indices[i][j];
			gridIndices[k++] = indices[i + 1][j + 1];
			gridIndices[k++] = indices[i + 1][j];

			velVertices[i][j][0] = Vertex((float)(j + 0.5) * cellSize, (float)(i + 0.5) * cellSize, 0.f, 0.f, 1.f, 1.f);
			velVertices[i][j][0] = Vertex((float)(j + 0.5) * cellSize, (float)(i + 0.5) * cellSize, 0.f, 0.f, 1.f, 1.f);
		}
	}

	updateGrid();
	simTime = 0;
}

/*
 * Unit square bilinear-interpolation (https://en.wikipedia.org/wiki/Bilinear_interpolation#On_the_unit_square).
 * You need this for interpolating the value in the cell (e.g., in backtracing advect).
 * Parameters:
 *	x, y:
 *	    position of the point whose value you want to interpolate
 *	v00:
 *	    value at (0, 0)
 *	v01:
 *	    value at (0, 1)
 *	v10:
 *	    value at (1, 0)
 *	v11:
 *	    value at (1, 1)
 *  return:
 *	interpolated value
 */
double bilinearInterpolate(double x, double y, double v00, double v01, double v10, double v11) {
	// TODO: Do a unit square bilinear interpolation for a point at (x, y).
	
    //a00 = v00;
	double a10 = v10 - v00;
	double a01 = v01 - v00;
	double a11 = v11 - v10 - v01 + v00;

	return v00 + a10*x + a01*y + a11*x*y;
}

/*
 * Set up the boundary.
 * Parameters:
 *	a:
 *	    2D array whose boundary need to be set
 *	flowType:
 *	    type of flow: horizontal flow (the interpolated quantity will goes to zero at the vertical boundaries), vertical flow (the interpolated quantity will goes to zero at the horizontal boundaries), other (e.g., temperature, density, which are not really flow... only continuity need to be guaranteed.)
 */
void setBnd(double a[rows + 2][cols + 2], FlowType flowType) {
	// TODO: Set up the boundary according to the flow type.
		


	for (int i = 1; i <= rows; i++) {
		a[i][0] = flowType == FlowType::horizontal ? -a[i][1] : a[i][1];
		a[i][cols + 1] = flowType == FlowType::horizontal ? -a[i][cols] : a[i][cols];

	}

	for (int j = 1; j <= cols; j++) {
		a[0][j] = flowType == FlowType::vertical ? -a[1][j] : a[1][j];
		a[rows + 1][j] = flowType == FlowType::vertical ? -a[rows][j] : a[rows][j];
	}

	a[0][0] = 0.5 * (a[0][1] + a[1][0]);
	a[rows + 1][0] = 0.5 * (a[rows + 1][1] + a[rows][0]);
	a[0][cols + 1] = 0.5 * (a[0][cols] + a[1][cols + 1]);
	a[rows + 1][cols + 1] = 0.5 * (a[rows + 1][cols] + a[rows][cols + 1]);

}

/*
 * Add source to the field.
 * Parameters:
 *	src:
 *	    2D array containing the source
 *	a:
 *	    target 2D array storing the quantity to modify
 */
void addSource(double src[rows + 2][cols + 2], double a[rows + 2][cols + 2]) {
	// TODO: Add the source from array *src* to the target array *a*: e.g., a = a + dt * src
	for (int i = 0; i < rows + 2; i++) {
		for (int j = 0; j < cols + 2; j++) {
			if (src[i][j] > 0) {
				a[i][j] = a[i][j] + dt * src[i][j];
			}
			else {

				a[i][j] = a[i][j] + dt * src[i][j];
			}

		}
	}


}

/*
 * Compute the diffusion part.
 * Parameters:
 *	a0:
 *	    2D array storing the quantities in current state
 *	a1:
 *	    2D array storing the quantities in next state
 *	nv:
 *	    diffusion rate (nv/kappa in the equations)
 *	flowType:
 *	    flow type
 */
void diffuse(double a0[rows + 2][cols + 2], double a1[rows + 2][cols + 2], double nv, FlowType flowType) {
	// TODO: diffusion
	// Compute the diffusion part and update array *a1*.
	// Use dt and dx for step size and cell size.
	// Do a implicit solve for stability.
	// Use a Gauss Seidel solve (or better, Conjugate Gradients).
	// Use *iterations* as the number of iterations.
	// Call setBnd to fix the boundary.

	double diff_rate = dt * nv * (rows * dx) * (cols * dx);

	for (int k = 0; k < iterations; k++) {
		for (int i = 1; i <= rows; i++) {
			for (int j = 1; j <= cols; j++) {
				a1[i][j] = (a0[i][j] + diff_rate * (a1[i - 1][j] + a1[i + 1][j] + a1[i][j - 1] + a1[i][j + 1])) / (1 + 4 * diff_rate);

			}
		}

		setBnd(a1, flowType);
	}



}

/*
 * Compute the advection part.
 * Parameters:
 *	a0:
 *	    2D array storing the quantities in current state
 *	a1:
 *	    2D array storing the quantities in next state
 *	vx, vy:
 *	    2D arrays storing the velocity field
 *	flowType:
 *	    flow type
 */
void advect(double a0[rows + 2][cols + 2], double a1[rows + 2][cols + 2], double vx[rows + 2][cols + 2], double vy[rows + 2][cols + 2], FlowType flowType) {
	// TODO: advection
	// Compute the advection part and update array *a1*.
	// Use dt and dx for step size and cell size.
	// Do a linear (or better, higher order or adaptive) backtrace for each center of the cells.
	// Compute the quantity in the previous state using the bilinear interpolation.
	// Call setBnd to fix the boundary.


	//double dt0 = dt / dx;
	double dt0x = dt * dx * cols;
	double dt0y = dt * dx * rows;
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			double x = j - dt0x * vx[i][j];
			double y = i - dt0y * vy[i][j];
			if (x < 0.5) x = 0.5; if (x > cols + 0.5) x = cols + 0.5; int j0 = (int)x; int j1 = j0 + 1;
			if (y < 0.5) y = 0.5; if (y > rows + 0.5) y = rows + 0.5; int i0 = (int)y; int i1 = i0 + 1;
			double s1 = x - j0; double s0 = 1 - s1; double t1 = y - i0; double t0 = 1 - t1;
			a1[i][j] = t0 * (s0 * a0[i0][j0] + s1 * a0[i0][j1]) + t1 * (s0 * a0[i1][j0] + s1 * a0[i1][j1]);
		}
	}

	setBnd(a1, flowType);

}


double s_p[rows + 2][cols + 2], s_div[rows + 2][cols + 2];

/*
 * Projection for the mass conservation.
 * Parameter:
 *	vx, vy:
 *	    the velocity field to be fixed
 */
void project(double vx[rows + 2][cols + 2], double vy[rows + 2][cols + 2]) {
	// TODO: projection
	// Do a Poisson Solve to get a divergence free velocity field.
	// Use a Gauss Seidel solve (or better, Conjugate Gradients).
	// Use *iterations* as the number of iterations.
	// Call setBnd to fix the boundary.

	double div[rows + 2][cols + 2];
	double p[rows + 2][cols + 2];

	memset(div, 0, sizeof(div));
	memset(p, 0, sizeof(p));


	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {

			div[i][j] = -0.5 * dx * (vy[i + 1][j] - vy[i - 1][j] + vx[i][j + 1] - vx[i][j - 1]);
			p[i][j] = 0;

		}
	}

	setBnd(div, FlowType::other);
	setBnd(p, FlowType::other);

	for (int k = 0; k < iterations; k++) {
		for (int i = 1; i <= rows; i++) {
			for (int j = 1; j <= cols; j++) {
				p[i][j] = (div[i][j] + p[i - 1][j] + p[i + 1][j] + p[i][j - 1] + p[i][j + 1]) / 4;
			}
		}
		setBnd(p, FlowType::other);
	}

	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			vy[i][j] -= 0.5 * (p[i + 1][j] - p[i - 1][j]) / dx;
			vx[i][j] -= 0.5 * (p[i][j + 1] - p[i][j - 1]) / dx;
		}
	}
	setBnd(vx, FlowType::horizontal);
	setBnd(vy, FlowType::vertical);

}

/*
 * Get the reference (average) temperature of the grid.
 * Parameters:
 *	tp:
 *	    2D array storing the temperatures
 */
double getReferenceTemperature(double tp[rows + 2][cols + 2]) {
	// TODO: Sum up array *tp* and compute the average.
	double tp_average = 0;
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			tp_average += tp[i][j];
		}
	}

 	tp_average = tp_average / (rows * cols);

	return tp_average;
}

/*
 * Apply the buoyancy force due to the temperature difference.
 * Parameters:
 *	a:
 *	    2D array storing the velocity
 *	tp:
 *	    2D array storing the temperature
 *	beta:
 *	    buoyancy force coefficient
 *	flowType:
 *	    flow type
 */
void applyTemperatureForce(double a[rows + 2][cols + 2], double tp[rows + 2][cols + 2], double beta, FlowType flowType) {
	// TODO: buoyancy forces
	// Apply the buoyancy force and update array *a*.
	// For more details, see Foster and Metaxas [1997] Equation 2.
	double ref = getReferenceTemperature(tp);

	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {

			a[i][j] = a[i][j] + dt * beta * (tp[i][j] - ref);
		}
	}
	
	setBnd(a, flowType);

}


void applyGravityForce(double a[rows + 2][cols + 2], double fg[rows + 2][cols + 2], double eg[rows + 2][cols + 2], double gravity, FlowType flowType) {

	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {

			a[i][j] = a[i][j] - dt * gravity * (fg[i][j] + eg[i][j]);
		}
	}

	setBnd(a, flowType);

}


void applyExpansionForce(double a[rows + 2][cols + 2], double tp[rows + 2][cols + 2], float reactionTemperature, double blueCoreCenter, double expansion, FlowType flowType) {
	
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			if (tp[i][j] > reactionTemperature) {
				a[i][j] = a[i][j] + dt * expansion * (j - blueCoreCenter);
			}
		}
	}

	setBnd(a, flowType);
}

void applyVorticityConfinement(double vx[rows + 2][cols + 2], double vy[rows + 2][cols + 2], double vorticityStrength) {

	double p[rows + 2][cols + 2];
	double p_diffx[rows + 2][cols + 2];
	double p_diffy[rows + 2][cols + 2];

	memset(p, 0, sizeof(p));
	memset(p_diffx, 0, sizeof(p_diffx));
	memset(p_diffy, 0, sizeof(p_diffy));

	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			p[i][j] = ((vx[i][j + 1] - vx[i][j - 1])/2 * dx) - ((vy[i + 1][j] - vy[i - 1][j]) / 2 * dx);
		}
	}

	//maybe abs the whole array instead
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			p_diffx[i][j] = ((abs(p[i][j + 1]) - abs(p[i][j - 1])) / 2 * dx);
			p_diffy[i][j] = ((abs(p[i + 1][j]) - abs(p[i + 1][j])) / 2 * dx);
			double norm = std::sqrt(p_diffx[i][j] * p_diffx[i][j] + p_diffy[i][j] * p_diffy[i][j]);
			
			if (norm > 0) {
				p_diffx[i][j] /= norm;
				p_diffy[i][j] /= norm;
			}
		}
	}

	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			vx[i][j] += dt * dx * vorticityStrength * p[i][j] * p_diffy[i][j];
			vy[i][j] -= dt * dx * vorticityStrength * p[i][j] * p_diffx[i][j];
		}
	}


	setBnd(vx, FlowType::horizontal);
	setBnd(vy, FlowType::vertical);

}


void combustion_decreaseFuel(double a[rows + 2][cols + 2], double tp[rows + 2][cols + 2], double fg[rows + 2][cols + 2], float reactionTemperature, double pctFuelBurn) {
	// TODO: Add the source from array *src* to the target array *a*: e.g., a = a + dt * src
	double burningRateParameter = -log(1 - pctFuelBurn);
	
	for (int i = 0; i < rows + 2; i++) {
		for (int j = 0; j < cols + 2; j++) {
			if (tp[i][j] > reactionTemperature) {
				a[i][j] = a[i][j] - dt * burningRateParameter * fg[i][j];
			}
		}
	}
}


void combustion_increaseExhaust(double a[rows + 2][cols + 2], double tp[rows + 2][cols + 2], double fg[rows + 2][cols + 2], float reactionTemperature, double pctFuelBurn, double stochiometricMixture) {
	// TODO: Add the source from array *src* to the target array *a*: e.g., a = a + dt * src
	double burningRateParameter = -log(1 - pctFuelBurn);
	double combustionParameter = burningRateParameter * stochiometricMixture;

	for (int i = 0; i < rows + 2; i++) {
		for (int j = 0; j < cols + 2; j++) {
			if (tp[i][j] > reactionTemperature) {
				a[i][j] = a[i][j] + dt * combustionParameter * fg[i][j] * (1 + 1 / burningRateParameter);
			}
		}
	}
}

void combustion_increaseTemperature(double a[rows + 2][cols + 2], double fg[rows + 2][cols + 2], float reactionTemperature, double pctFuelBurn, double stochiometricMixture, double combustionHeat) {
	// TODO: Add the source from array *src* to the target array *a*: e.g., a = a + dt * src
	double burningRateParameter = -log(1 - pctFuelBurn);
	double combustionParameter = burningRateParameter * stochiometricMixture;

	for (int i = 0; i < rows + 2; i++) {
		for (int j = 0; j < cols + 2; j++) {
			if (a[i][j] > reactionTemperature) {
				a[i][j] = a[i][j] + dt * combustionParameter * fg[i][j] * combustionHeat;
			}
		}
	}
}


double planckFormula(double tp, double wv) {
	double numerator = 2 * pi * planck_cst * light_spd * light_spd;
	double denominator = pow(wv, 5) * (exp((planck_cst * light_spd) / (wv * tp * kBoltzmann)) - 1);
	return numerator / denominator;
}

double mapIntensities(double L, double Laverage) {
	return 1 * (1 - exp(-L / Laverage));
}


void applyBlackBodyRadiation(double tp[rows + 2][cols + 2]) {
	for (size_t i = 0; i <= rows; ++i)
		for (size_t j = 0; j <= cols; ++j) {

			float v00 = tp[i][j];
			float v01 = tp[i + 1][j];
			float v10 = tp[i][j + 1];
			float v11 = tp[i + 1][j + 1];
			double t = bilinearInterpolate(0.5, 0.5, v00, v01, v10, v11);
			gridVertices[i][j].r = mapIntensities(planckFormula(t, wavelenghtRed), mapIntensity);
			gridVertices[i][j].g = mapIntensities(planckFormula(t, wavelenghtGreen), mapIntensity);
			gridVertices[i][j].b = mapIntensities(planckFormula(t, wavelenghtBlue), mapIntensity);

		}

}



void applyBlackBodyRadiationSmoke(double tp[rows + 2][cols + 2], double eg[rows + 2][cols + 2]) {
	for (size_t i = 0; i <= rows; ++i)
		for (size_t j = 0; j <= cols; ++j) {
			if (eg[i][j] > 0) {
				float v00 = tp[i][j];
				float v01 = tp[i + 1][j];
				float v10 = tp[i][j + 1];
				float v11 = tp[i + 1][j + 1];
				double t = bilinearInterpolate(0.5, 0.5, v00, v01, v10, v11);
				gridVertices[i][j].r = eg[i][j] * mapIntensities(planckFormula(t, wavelenghtRed), mapIntensity);
				gridVertices[i][j].g = eg[i][j] * mapIntensities(planckFormula(t, wavelenghtGreen), mapIntensity);
				gridVertices[i][j].b = eg[i][j] * mapIntensities(planckFormula(t, wavelenghtBlue), mapIntensity);
			}

		}

}



void dissipate(double a[rows + 2][cols + 2], double dissipation) {
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			if (a[i][j] > 0) {
				a[i][j] = a[i][j] - dissipation * dt * a[i][j];
			}
		}
	}
}







void containSource(double src[rows + 2][cols + 2], double a[rows + 2][cols + 2]) {
	for (int i = 1; i <= rows; i++) {
		for (int j = 1; j <= cols; j++) {
			if (src[i][j] != 0) {
				a[i][j] = src[i][j];
			}
		}
	}
}


void draw_race(int start) {	
	
	double vel = (double) rand() / RAND_MAX;

	vel = vel <= 0.3 ? 0.3 : vel;

	for (int i = 10; i <= 32; i++) {
		cur_vx[start][i] = vel;
	}

	for (int j = start-5; j <= start; j++) {
		cur_vy[j][32] = -vel;
	}

	for (int i = 32; i <= 42; i++) {
		cur_vx[start-5][i] = vel;
	}

	for (int j = start-5; j <= start + 7; j++) {
		cur_vy[j][42] = vel;
	}

	for (int i = 32; i <= 42; i++) {
		cur_vx[start + 7][i] = -vel;
	}

	for (int j = start + 7; j <= start + 11; j++) {
		cur_vy[j][32] = vel;
	}

	for (int i = 32; i <= 66; i++) {
		cur_vx[start + 11][i] = vel;
	}

	for (int j = start -4; j <= start + 11; j++) {
		cur_vy[j][66] = -vel;
	}

	for (int i = 66; i <= 90; i++) {
		cur_vx[start - 4][i] = vel;
	}

}



/*
 * One stimulation step.
 */

void snake_init() {
	snake_flag = true;

	kappa = 0;

	int start1 = 7;
	tpSrc[start1][10] = 5000000;
	draw_race(start1);
	tpSrc[start1-5][90] = -100;
	tpSrc[start1-3][90] = -100;
	tpSrc[start1-4][91] = -100;

	int start2 = 27;
	tpSrc[start2][10] = 5000000;
	draw_race(start2);
	tpSrc[start2 - 5][90] = -100;
	tpSrc[start2 - 3][90] = -100;
	tpSrc[start2 - 4][91] = -100;

	int start3 = 47;
	tpSrc[start3][10] = 5000000;
	draw_race(start3);
	tpSrc[start3 - 5][90] = -100;
	tpSrc[start3 - 3][90] = -100;
	tpSrc[start3 - 4][91] = -100;

	int start4 = 67;
	tpSrc[start4][10] = 5000000;
	draw_race(start4);
	tpSrc[start4 - 5][90] = -100;
	tpSrc[start4 - 3][90] = -100;
	tpSrc[start4 - 4][91] = -100;

}



void step() {

	if (snake_flag){
		addSource(tpSrc, next_tp);

		std::swap(cur_tp, next_tp);
		diffuse(cur_tp, next_tp, kappa, FlowType::other);
		std::swap(cur_tp, next_tp);
		advect(cur_tp, next_tp, cur_vx, cur_vy, FlowType::other);

	}

	else {
		// TODO: step the simulation
		// Compute and update the velocity and temperature (in cur_vx, cur_vy, cur_tp) in the next state based on the current state.
		// You need to apply source, diffuse, advect forces for temperature and velocity fields.
		// For velocity field, you need to do the projection to get a divergence free field.
		// You also need to apply the buoyancy force to the velocity field.
		// Don't forget to swap pointers (e.g., cur_vx and next_vx, etc.)!
		// Have a look at GDC03.
		// Change the paramters (dt, dx, etc.) in the setting panel to check whether your solver is stable!



		//===VEL STEP===
		addSource(vxSrc, cur_vx);
		addSource(vySrc, cur_vy);

		std::swap(cur_vx, next_vx);
		std::swap(cur_vy, next_vy);

		diffuse(cur_vx, next_vx, nv, FlowType::horizontal);
		diffuse(cur_vy, next_vy, nv, FlowType::vertical);

		applyVorticityConfinement(next_vx, next_vy, vorticityStrength);
		applyTemperatureForce(next_vy, cur_tp, buoyancy, FlowType::vertical); //buoyancy
		applyGravityForce(next_vy, cur_fg, cur_eg, gravity, FlowType::vertical);
		applyExpansionForce(next_vx, cur_tp, reactionTemperature, blueCoreCenter, expansion, FlowType::horizontal);


		project(next_vx, next_vy);

		std::swap(cur_vx, next_vx);
		std::swap(cur_vy, next_vy);

		advect(cur_vx, next_vx, cur_vx, cur_vy, FlowType::horizontal);
		advect(cur_vy, next_vy, cur_vx, cur_vy, FlowType::vertical);

		project(next_vx, next_vy);




		//===DENSITY STEP===
		dissipate(next_tp, dissipation);
		dissipate(next_fg, dissipation_fuel);
		dissipate(next_eg, dissipation_exhaust);



		addSource(tpSrc, next_tp);
		addSource(fgSrc, next_fg);

		//COMPUTE COMBUSTION
		combustion_increaseExhaust(next_eg, cur_tp, cur_fg, reactionTemperature, pctFuelBurn, stochiometricMixture);
		combustion_increaseTemperature(next_tp, cur_fg, reactionTemperature, pctFuelBurn, stochiometricMixture, combustionHeat);
		combustion_decreaseFuel(next_fg, cur_tp, cur_fg, reactionTemperature, pctFuelBurn);
		

		std::swap(cur_tp, next_tp);
		diffuse(cur_tp, next_tp, kappa, FlowType::other);
		std::swap(cur_tp, next_tp);
		advect(cur_tp, next_tp, next_vx, next_vy, FlowType::other);


		//===FUEL STEP===


		std::swap(cur_fg, next_fg);
		diffuse(cur_fg, next_fg, kappa_fuel, FlowType::other);
		std::swap(cur_fg, next_fg);
		advect(cur_fg, next_fg, next_vx, next_vy, FlowType::other);


		//===EXHAUST STEP===


		std::swap(cur_eg, next_eg);
		diffuse(cur_eg, next_eg, kappa_exhaust, FlowType::other);
		std::swap(cur_eg, next_eg);
		advect(cur_eg, next_eg, next_vx, next_vy, FlowType::other);






	}

	// Please DO NOT change the following
	simTime += dt;
	//updateGrid();
	updateFilament();
	memset(vxSrc, 0, sizeof(vxSrc));
	memset(vySrc, 0, sizeof(vySrc));
}

void updateGrid() {
	if (vizType == 0) {
		applyBlackBodyRadiation(cur_tp);
	}
	else if(vizType == 1) {
		applyBlackBodyRadiation(cur_fg);
	}	
	else{
		//ONLY CONSIDER cur_tp where there is some exhaust for this one
		applyBlackBodyRadiationSmoke(cur_tp, cur_eg);
	}

	for (size_t i = 1; i <= rows; ++i)
		for (size_t j = 1; j <= cols; ++j) {
			velVertices[i - 1][j - 1][1].x = velVertices[i - 1][j - 1][0].x + cur_vx[i][j] * velScale;
			velVertices[i - 1][j - 1][1].y = velVertices[i - 1][j - 1][0].y + cur_vy[i][j] * velScale;
		}
}

double dist(const Vertex &a, const Vertex &b) {
    double x = a.x - b.x;
    double y = a.y - b.y;
    return std::sqrt(x*x + y*y);
}

void addFilament(double x, double y) {
	filaments.push_back(make_shared<vector<Vertex>>());
	for (size_t i = 0; i <= rows; ++i) {
		filaments.back()->push_back(Vertex(x, (float)i * cellSize, 1.0f, 1.0f, 1.0f, 1.0f));
	}
	ages.push_back(0);

	filaments.push_back(make_shared<vector<Vertex>>());
	for (size_t i = 0; i <= cols; ++i) {
		filaments.back()->push_back(Vertex((float)i * cellSize, y, 1.0f, 1.0f, 1.0f, 1.0f));
	}
	ages.push_back(0);
}

void updateFilament() {
	while (!ages.empty() && ages.front() > maxAge) {
		ages.pop_front();
		filaments.pop_front();
	}
	for (size_t i = 0; i < filaments.size(); ++i) {
		ages[i] += dt;
		shared_ptr<vector<Vertex>> tmp = make_shared<vector<Vertex>>();
		tmp->push_back(filaments[i]->front());
		for (size_t j = 1; j < filaments[i]->size(); ++j) {
			const Vertex &a = tmp->back();
			const Vertex &b = filaments[i]->at(j);
			if (dist(a, b) > refinementThreshold) {
				tmp->push_back(Vertex((a.x + b.x)/2, (a.y + b.y)/2, b.r, b.g, b.b, b.a));
			}
			tmp->push_back(b);
		}
		filaments[i] = tmp;
		for (Vertex &v: *filaments[i]) {
			double x = v.x/cellSize + 0.5;
			double y = v.y/cellSize + 0.5;
			size_t j0 = (size_t)x;
			size_t i0 = (size_t)y;
			size_t i1 = i0 + 1, j1 = j0 + 1;
			x -= j0;
			y -= i0;
			double v00, v01, v10, v11;
			double vx, vy;
			v00 = cur_vx[i0][j0];
			v01 = cur_vx[i1][j0];
			v10 = cur_vx[i0][j1];
			v11 = cur_vx[i1][j1];
			vx = bilinearInterpolate(x, y, v00, v01, v10, v11);
			v00 = cur_vy[i0][j0];
			v01 = cur_vy[i1][j0];
			v10 = cur_vy[i0][j1];
			v11 = cur_vy[i1][j1];
			vy = bilinearInterpolate(x, y, v00, v01, v10, v11);
			v.x += (float)vx * dt * cellSize;
			v.y += (float)vy * dt * cellSize;
			v.x = std::max(0.f, v.x);
			v.x = std::min((float)frameWidth, v.x);
			v.y = std::max(0.f, v.y);
			v.y = std::min((float)frameHeight, v.y);
			if (ages[i] > maxAge / 2)
			v.a = 2 - (float)ages[i] / (maxAge / 2);
		}
	}
}
