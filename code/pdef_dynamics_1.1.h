/*
 
  Code for the article 'Animal diversity and ecosystem functioning in dynamic food webs'
   
  Copyright (C) 2016 Christian Guill & Florian D. Schneider
  
  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
  
   
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
