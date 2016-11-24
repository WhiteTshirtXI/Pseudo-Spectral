#ifndef _PAEAMETERS_
#define _PAEAMETERS_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <algorithm>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//#define GHOST_BOUNDARIES

#define algorithm FFTW_EXHAUSTIVE
//#define algorithm FFTW_MEASURE
//#define algorithm FFTW_PATIENT


using namespace std;

#define index(x,y) (x)*Ny+y

typedef double Real;

enum ComputeMode {tr, ga};
ComputeMode compute_mode;

const int Nx = 128; // This number must be in the form of 2k (because we are using image method)
const int Ny = Nx/4;
const int N = Nx*Ny;
#ifdef GHOST_BOUNDARIES
const int Nx_box = Nx / 2 + 1;
const int Ny_box = Ny / 2 + 1;
#else
const int Nx_box = Nx;// Nx / 2 + 1;
const int Ny_box = Ny;//Ny / 2 + 1;
#endif
Real Lx_box = 64; // Box x is from -Lx to Lx
Real Ly_box = (Lx_box*Ny) / Nx; // Box y is from -Ly to Ly
Real Lx = Nx*Lx_box / Nx_box;
Real Ly = (Ny*Ly_box) / Ny_box;
const int Nkx = Nx;
const int Nky = Ny/2 + 1;
const int Nk = Nkx*Nky;
Real kx[Nkx];
Real ky[Nky];

const int aao = 2;// anti aliasing order
const int kcy = 2*Nky/(aao+1);
const int kcx = 2*(Nx/2+1)/(aao+1);
const int default_interval = 16; // interval is the number of step of simulation in unit time.
Real dt = 1.0/(default_interval);
const Real eta = 0.5;
const Real epsilone = 0.01;//dt*dt*dt*dt;

const double time_digits = 100000.0;

const gsl_rng_type * gsl_T;
gsl_rng * gsl_r;

Real initial_rho = 8;
Real T = 1;
Real v = 1;
Real g = 0.125;
Real alpha = 0.0;
Real E = v/2;
const int layer = 0;
long double lapsed_time = 0;

Real K = 0.125;
Real D = K;

#endif

