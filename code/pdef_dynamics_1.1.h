/*
 
 MIT License
 
 Copyright (C) 2016 Christian Guill & Florian D. Schneider
 
 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:
 
 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.
 
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <stdlib.h>	
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>				// random number generator
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_blas.h>				// linear algebra routines
#include <gsl/gsl_sort_vector.h>		// vector operations
#include <gsl/gsl_odeiv.h>              // ODE solver 

#include <cvode/cvode.h>             /* prototypes for CVODE fcts., consts. */
#include <nvector/nvector_serial.h>  /* serial N_Vector types, fcts., macros */
#include <cvode/cvode_dense.h>       /* prototype for CVDense */


double *web_calc(gsl_rng *r,FILE *data1);
static void dynamics(double B[], double Bdot[], void *params);
static void pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mass);
static void output(gsl_matrix *Ap,gsl_vector *mas,double B[],double iniB[],double meanB[],double VarCoeffB[],void *params,FILE *data1);
static void Average_TL(gsl_matrix *Ap, double TLav[], double TLmin[]);
static void Min_TL(gsl_matrix *Ap, double TLmin[]);
static void Eff_TL(gsl_matrix *Ap, double meanB[], void *params, double TLeff[], double TLmin[]);
double functional_diversity(gsl_vector *mas);


double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff ,gsl_rng *r);
static void mean_body_mass(double meanB[], gsl_vector *mass, double mean_body_masses[]);

static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, double params[]);
static void initialise_biomass_densities(double B[], double iniB[], double params[], gsl_rng *r);
static void solve_ode(double B[], double meanB[], double CV[], double params[]);
static void Extinction(double B[], int N);
static int Cvode_rates(realtype t, N_Vector y, N_Vector ydot, void *params);

static void net_rates(gsl_vector *TBvec, gsl_vector *Ivec, gsl_vector *Dvec, gsl_vector *Xvec, gsl_vector *N_in, gsl_vector *N_out, gsl_vector *N_up, void *params);
static void flows(double meanB[], void *params, double energyflows[]);
static void mass_abundance(double meanB[], gsl_vector *mas, double fitresults[]);

static void resolve_input(int argc, char *input[]);
static void write_header(FILE *data1);
static void Prepare_timeseries_file(FILE *timeseries);
static void Write_timeseries_to_file(double B[], double t, FILE *timeseries);

static void show_matrix(gsl_matrix *A, int N);
