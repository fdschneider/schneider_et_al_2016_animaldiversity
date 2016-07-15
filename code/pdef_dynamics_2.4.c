/*
Project: Predator Diversity and Ecosystem Functioning
------
Population Dynamics Code

version 2.4.1
last edit: 14.07.2016 by CG

authors: Christian Guill & Florian D. Schneider

comments:  * program has to be invoked with a parameter to seed the RNG (two-digit number, i.e., 01)
           * updated version following the completion of the revision process; simulations were using version 2.4 at commit #74cbe2dd)
*/

#include "pdef_dynamics_1.1.h"

//  basic simulation variables:
int Web_Id=0;
int S;                        // total number of species
int S_c=30;                   // number consumers (animals, predators)
int S_b=30;                   // number basal species (plants)
int N=2;                      // number of nutrients

double tend=150000;           // time at which simulation runs stop
double teval=10000;           // length of time interval from which to determine mean(B) and CV(B)
double Delta_t = 0.1;         // step size in ODE integration (internal step size can be less or larger) 

const double eps_rel= 1e-10;  // relative error tolerance for ODE solver
const double eps_abs = 1e-10; // absolute error tolerance for ODE solver
const double EXTINCT = 1e-6;  // extinction threshold for biomass densities
int seed=450000;              // seed for random number generator (modified by command-line input)
int Imax=40;                  // number of replicates per diversity level
int DivLev=100;               // number of diversity levels
int DivStep=1;                // step size of diversity levels
char filename[256]="pdef_2.4_H50_";          // name of output file
int job_ID;                                  // will be read in as a command-line parameter
char timeseries_file[256]="timeseries.out";  // name of output file for the time series (will be modified later)

//  parameters that determine the network topology
double zeta_c=6;              // log_10 of the mean of the consumer (animal) body masses
double sigma_c=3;             // width of the distribution of consumer (animal) body masses
double cutoff_c=1e5;          // half relative cutoff for distribution of consumer body masses
double zeta_b=5;              // log_10 of the mean of the basal (plant) body masses
double sigma_b=3;             // width of the distribution of basal (plant) body masses
double cutoff_b=1e5;          // half relative cutoff for distribution of basal body masses

double m_p_min = 0;           // minimla and maximal log10 body masses in case of uniform distributions
double m_a_min = 2;
double m_p_max = 6;
double m_a_max = 12;

double cutoff=0.01;           // cutoff of the Ricker curve for setting a link between predator and prey
double R_opt=100;             // optimal predator-prey body-mass ratio
double g=2;                   // width of the Ricker curve (higher g-value -> narrower curve)

double f_herbiv = 0.50;       // fraction of species that are strict herbivores
double f_pred = 0.0;          // fraction of species that are strict predators

//  parameters of the functional response
double a_0=50;                // scaling factor for the attack rate
double a_c;                   // exponent for predator body-mass scaling of the attack rate (value is drawn from distribution)
double mean_a_c = 0.47;       // (mean for the above distribution)
double sigma_a_c = 0.04;      // (standard deviation for the above distribution)
double a_p;                   // exponent for prey body-mass scaling of the attack rate (value is drawn from distribution)
double mean_a_p = 0.15;       // (mean for the above distribution)
double sigma_a_p = 0.03;      // (standard deviation for the above distribution)
double a_plant = 20;          // constant attack rate coefficient for plants

double h_0=0.4;               // scaling factor for the handling time
double h_c;                   // exponent for predator body-mass scaling of the handling time (value is drawn from distribution)
double mean_h_c = -0.48;      // (mean for the above distribution)
double sigma_h_c = 0.03;      // (standard deviation for the above distribution)
double h_p;                   // exponent for prey body-mass scaling of the handling time (value is drawn from distribution)
double mean_h_p = -0.66;      // (mean for the above distribution)
double sigma_h_p = 0.02;      // (standard deviation for the above distribution)

double hill;                  // Hill coefficient (value is drawn from distribution
double mu_hill=1.5;           // (mean for the above distribution)
double sigma_hill=0.2;        // (standard deviation for the above distribution)
double C_interference;        // interference competition (value is drawn from distribution)
double mu_Cint=0.8;           // (mean for the above distribution)
double sigma_Cint=0.2;        // (standard deviation for the above distribution)

double x_resp_a = 0.314;      // intercept of animal respiration rates
double x_resp_p = 0.138;      // intercept of producer respiration rates
double e_p=0.45;              // assimilation efficiency for plant resources
double e_a=0.85;              // assimilation efficiency for animal resources

double B_b=10;                // intercept initial basal biomass densities
double B_c=10;                // intercept initial consumer biomass densities

//  parameters of the nutrient model
double C1=1;                  // content of first nutrient in plants
double C2=0.5;                // content of second nutrient in plants
double D_N=0.25;              // nutrient turnover rate
double K_min=0.1;             // minimal nutrient uptake half saturation density
double K_rel=0.1;             // width of interval for nutrient uptake half saturation densities
double mu_S=10;               // mean nutrient supply concentration
double sigma_S=2;             // standard deviation of nutrient supply concentration

//  default parameters for loops and flags
int GF=0;                     // global flag, for testing purposes
double cv_flag;               // indicates potentially pathological coefficients of variation
int TIMESERIES = 0;           // if true, for the last network a file with the time series will be created  
int VAR_COEFF = 1;            // if true, the coefficient of variation will be calculated; if false, the simulations run a bit faster
int UNIFORM_MASSES = 1;       // 0: body masses drawn from log-normal distribution; 1: log-body masses drawn from uniform distribution

//  ********** Main function **********
int main(int argc, char* argv[]) {
  FILE *data1;

  int i,j;

  resolve_input(argc,argv);                 // read in the command-line parameter, modify rng seed and output-file name
  Web_Id=DivLev*Imax*(job_ID-1);            // increase Web_Id

  gsl_rng_default_seed = seed;              // seed the rng
  gsl_rng *r=gsl_rng_alloc(gsl_rng_default);

  write_header(data1);                      // write a header for the output file

  for(j=10; j<=DivLev; j++) {
    S_c=j*DivStep;
    S=S_c+S_b;

    for(i=0;i<Imax;i++) {
      printf("%d  %d\n",j,i);
      Web_Id++;

      data1=fopen(filename,"a");            //re-opens output file to amend it
      web_calc(r,data1);
      fclose(data1);                        // forces programm to write data to file
    }
  }

  gsl_rng_free(r);

  return(0);
}

//  ********** generate a network, simulate population dynamics, and evaluate the final network **********
double *web_calc(gsl_rng *r,FILE *data1) {
  gsl_matrix *Ap=gsl_matrix_calloc(S,S);                   // adjacency matrix
  gsl_vector *mass=gsl_vector_calloc(S);                   // mean body masses of the species

  double *B=(double *) calloc(S+N,sizeof(double));         // biomass densities
  double *iniB=(double *) calloc(S+N,sizeof(double));      // initial biomass denisties
  double *meanB=(double *) calloc(S+N,sizeof(double));     // mean biomass densities 
  double *VarCoeffB=(double *) calloc(S+N,sizeof(double)); // coefficients of variation
  double *params=(double *) calloc(2*S*S+6*S+S_b*N+N+S_b+N,sizeof(double));  // parameters for population dynamics

  pdef_structure(r,Ap,mass);                               // generate network structure
  set_parameters(r,Ap,mass,params);                        // set parameters for the simulation run
  initialise_biomass_densities(B, iniB, params, r);        // set initial values of biomass densities
  solve_ode(B,meanB,VarCoeffB,params);                     // apply solver for integration of the ODE
  output(Ap,mass,B,iniB,meanB,VarCoeffB,params,data1);     // evaluate network and generate output

//  *********** Free memories ***********
  gsl_matrix_free(Ap);
  gsl_vector_free(mass);
  free(params);
  free(meanB);
  free(VarCoeffB);

  return 0;
}

//  ********** generate the body-mass based network structure **********
static void pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mass) {
  int i,j,flag=0;
  double Wkeit,sigma_i,zeta_act;
  double temp1,temp2,m_opt,a_max,R;
  double m_crit;

  while(flag==0) {
    flag=1;
    gsl_vector_set_all(mass,0);
    gsl_matrix_set_zero(Ap);

//  ********* determine body masses *********
    if(UNIFORM_MASSES == 0) {
      zeta_act=zeta_c*log(10);                                          // calculate log(mean) from log10(mean)
      for(i=0;i<S_c;i++) {
        temp1=gsl_ran_lognormal(r,zeta_act,sigma_c);  
        if(temp1>exp(zeta_act)/cutoff_c&&temp1<exp(zeta_act)*cutoff_c)  // check whether mass is within certain boundaries
          gsl_vector_set(mass,S_b+i,temp1);    
        else
          i--;
      }
      gsl_sort_vector(mass);

      gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
      gsl_vector *mass_b=&mass_b_vec.vector;

      zeta_act=zeta_b*log(10);                                          // calculate log(mean) from log10(mean)
      for(i=0;i<S_b;i++) {
        temp1=gsl_ran_lognormal(r,zeta_act,sigma_b);  
        if(temp1>exp(zeta_act)/cutoff_b&&temp1<exp(zeta_act)*cutoff_b)  // check whether mass is within certain boundaries
          gsl_vector_set(mass_b,i,temp1);              
        else
          i--;
      }
      gsl_sort_vector(mass_b);
    }
    else {
      for(i = 0; i< S_c; i++)
        gsl_vector_set(mass,S_b+i,pow(10.,gsl_ran_flat(r,m_a_min,m_a_max)));
      gsl_sort_vector(mass);

      gsl_vector_view mass_b_vec=gsl_vector_subvector(mass,0,S_b);
      gsl_vector *mass_b=&mass_b_vec.vector;
      
      for(i = 0; i< S_b; i++)
        gsl_vector_set(mass_b,i,pow(10.,gsl_ran_flat(r,m_p_min,m_p_max)));
      gsl_sort_vector(mass_b);
    }

//  ********** fill the adjacency matrix with Ricker attack rates *************
    for(i=0;i<S_c;i++) {
      temp1=gsl_vector_get(mass,S_b+i);
      a_max=pow(gsl_vector_get(mass,S_b+i)/R_opt,0.25);

      for(j=0;j<S;j++) {
        temp2=gsl_vector_get(mass,j);
        R = temp1/temp2;
        if(pow((R / R_opt) * exp(1 - (R / R_opt)), g) >= cutoff)
          gsl_matrix_set(Ap,S_b+i,j,1);
      }

      if(gsl_rng_uniform(r) < (f_herbiv + f_pred)) {                    // combined probability to mess with this species' links
        if(gsl_rng_uniform(r) < f_herbiv/(f_herbiv+f_pred)) {           // either make it a strict herbivore...
          if(UNIFORM_MASSES == 0)
            m_crit = pow(10,zeta_b) * cutoff_b * R_opt;
          else
            m_crit = pow(10,m_p_max) * R_opt;
            
          if(gsl_vector_get(mass,S_b+i) < m_crit) {    
            for(j = S_b; j < S; j++)
              gsl_matrix_set(Ap,S_b+i,j,0);                             // and remove all links from non-plant resources
          }
        }
        else {                                                          // ... or make it a strict carnivore 
          for(j = 0; j < S_b; j++)
            gsl_matrix_set(Ap,S_b+i,j,0);                               // and remove all links from plant resources
        }
      }

      gsl_vector_view tempp=gsl_matrix_row(Ap,S_b+i);                   // reject networks with consumers or predators without prey
      flag=flag*(1-gsl_vector_isnull(&tempp.vector));
    }
    
    for(i=0; i<S_b; i++) {
      gsl_vector_view tempp = gsl_matrix_column(Ap,i);                  // reject networks with uncontrolled basal species
      flag = flag*(1-gsl_vector_isnull(&tempp.vector));
    }
  }

  return;
}

//  ********** write all parameters required for the dynamics to the array 'params' **********
static void set_parameters(gsl_rng *r, gsl_matrix *Ap, gsl_vector *mass, double params[]) {
  int i,j;
  double temp1,temp2,temp3,temp4,R;

  gsl_matrix *A=gsl_matrix_calloc(S,S);                                 // attack rates
  gsl_matrix *H=gsl_matrix_calloc(S,S);                                 // handling times

  gsl_matrix_memcpy(A,Ap);

  hill = get_random_parameter(mu_hill,sigma_hill,1,2,r);
  C_interference = get_random_parameter(mu_Cint,sigma_Cint,mu_Cint - 3*sigma_Cint, mu_Cint + 3*sigma_Cint,r);
  
  a_c = get_random_parameter(mean_a_c, sigma_a_c, mean_a_c - 3*sigma_a_c, mean_a_c + 3*sigma_a_c,r);
  a_p = get_random_parameter(mean_a_p, sigma_a_p, mean_a_p - 3*sigma_a_p, mean_a_p + 3*sigma_a_p,r);

  h_c = get_random_parameter(mean_h_c, sigma_h_c, mean_h_c - 3*sigma_h_c, mean_h_c + 3*sigma_h_c,r);
  h_p = get_random_parameter(mean_h_p, sigma_h_p, mean_h_p - 3*sigma_h_p, mean_h_p + 3*sigma_h_p,r);
  

  for(i=0;i<S;i++) {
    if(i >= S_b) {                                                      // the following lines are only for non-basal species
      gsl_vector_view tempp=gsl_matrix_row(A,i);  
      temp1=gsl_blas_dasum(&tempp.vector);                              // temp1 stores the number of prey species of predator i
      gsl_vector_scale(&tempp.vector,1/temp1);                          // reduce attack rates for generalists 

      temp1 = pow(gsl_vector_get(mass,i),a_c);

      for(j=0;j<S;j++) {
        if(j < S_b)
          temp2 = a_plant;
        else
          temp2 = pow(gsl_vector_get(mass,j),a_p);
        R = gsl_vector_get(mass,i) / gsl_vector_get(mass,j);            // predator-prey body-mass ratio for Ricker curve
        temp3 = temp2 * temp1 * pow((R / R_opt) * exp(1 - (R / R_opt)),g);
        gsl_matrix_set(A,i,j,a_0*temp3*gsl_matrix_get(A,i,j));          // attack rates
        
        temp4 = h_0 * pow(gsl_vector_get(mass,i),h_c) * pow(gsl_vector_get(mass,j),h_p);
        gsl_matrix_set(H,i,j,temp4);                                    // handling times
      }
    }

    for(j=0;j<S;j++) {
      *(params+i*S+j)=gsl_matrix_get(A,i,j);
      *(params+S*S+i*S+j)=gsl_matrix_get(H,i,j);
    }

    if(i < S_b) {  
      *(params+2*S*S+i)=x_resp_p*pow(gsl_vector_get(mass,i),-0.25);     // plant respiration rates
      *(params+2*S*S+2*S+i)=1;                                          // identifyer for basal species
      *(params+2*S*S+5*S+i) = e_p;                                      // assimilation efficiency for plant resources
    }
    else {
      *(params+2*S*S+i)=x_resp_a*pow(gsl_vector_get(mass,i),-0.25);     // animal respiration rates
      *(params+2*S*S+5*S+i) = e_a;                                      // assimilation efficiency for animal resources
    }

    *(params+2*S*S+S+i)=C_interference;                                 // interference competition
    *(params+2*S*S+3*S+i)=hill;                                         // Hill coefficient
    *(params+2*S*S+4*S+i)=gsl_vector_get(mass,i);                       // body masses    
  }

  for(i=0;i<S_b;i++) {                                                  // the following lines are only for basal species
    for(j=0;j<N;j++)
      *(params+2*S*S+6*S+i*N+j)=K_min+K_rel*gsl_rng_uniform(r);         // nutrient uptake half saturation densities
  }
  
  for(i=0;i<N;i++) {
    temp1=gsl_ran_gaussian(r,sigma_S)+mu_S;
    if(temp1<=0)
      i--;
    else
      *(params+2*S*S+6*S+S_b*N+i)=temp1;                                // nutrient supply concentrations
  }
  for(i=0;i<S_b;i++)
    *(params+2*S*S+6*S+S_b*N+N+i)=pow(gsl_vector_get(mass,i),-0.25);    // max. nutrient uptake rates
    
  *(params+2*S*S+6*S+S_b*N+N+S_b) = C1;
  *(params+2*S*S+6*S+S_b*N+N+S_b+1) = C2;

  gsl_matrix_free(A);
  gsl_matrix_free(H);
  
  return;
}

// ********** get a random number from a Gaussian distribution with specified parameters and within specified limits **********
double get_random_parameter(double mean, double sigma, double low_cutoff, double high_cutoff ,gsl_rng *r) {
  double par = low_cutoff - 1;                                          // initialise parameter with a value outside the desired range
  
  while(par < low_cutoff || par > high_cutoff)
    par = gsl_ran_gaussian(r,sigma) + mean;

  return par;
}
  
//  ********** set the initial biomass densities and save them in an extra array **********
static void initialise_biomass_densities(double B[], double iniB[], double params[], gsl_rng *r) {
  int i;
  double temp1;
    
  for(i=0;i<S_b;i++)
    B[i]=B_b*gsl_rng_uniform_pos(r);                                    // basal species

  for(i=S_b;i<S;i++)
    B[i]=B_c*gsl_rng_uniform_pos(r);                                    // consumers and predators

  for(i=0;i<N;i++) {
    temp1=(double) *(params+2*S*S+5*S+2*S_b*N+i);
    B[S+i]=0.5*temp1+0.5*temp1*gsl_rng_uniform_pos(r);                  // nutrients
  }
  for(i=0;i<S+N;i++)
    iniB[i]=B[i];                                                       // saving initial biomass values
    
  return;
}

//  ********** read the input parameter and modify program parameters with it **********
static void resolve_input(int argc, char *input[]) {
  if(argc != 2) {
    printf("Please start job with an integer parameter\n");             // print a warning if no input parameter is given
    tend = 1e-7;
    teval = 1e-7;
  }

  job_ID = atoi( input[1]);                                             // read input parameter as integer
  const char *str = input[1];                                           // read it again as a character array

  seed += job_ID;                                                       // modify rng seed with the input parameter
  strcat(filename,str);                                                 // attach the parameter to the name of the output file
  
  return;
}

//  ********** prepare the output file and write the column headers (used in R) to the file **********
static void write_header(FILE *data1) {
  int i;
  data1=fopen(filename,"w");                                            // creates empty file with the name given in the header of the code

  fprintf(data1,"#job_id=%d,	seed=%d\n",job_ID,seed);
  fprintf(data1,"WebID	S_b	S_c	hill	pred_int	a_c	a_p	h_c	h_p	P	P_b	P_c	P_n	");
  fprintf(data1,"meanB	meanB_b	meanB_c	meanB_n	CV	CV_b	CV_c	CV_n	");
  fprintf(data1,"FD_ini	FD_fin	inidens_b	inidens_c	findens_b	findens_c	");
  for(i=0;i<N;i++)
    fprintf(data1,"nut_S%d	",i+1);
  fprintf(data1,"nut_in	nut_bas	bas_cons	cons_igp	bas_resp	cons_resp	");
  fprintf(data1,"meanTLmin_ini	meanTLmin_fin	meanTLav_ini	meanTLav_fin	meanTLeff	");
  fprintf(data1,"inilogmas_b	inilogmas_c	inivul_b	inivul_c	inigen_c	");	
  fprintf(data1,"finlogmas_b	finlogmas_c	finvul_b	finvul_c	fingen_c	");
  fprintf(data1,"iniSDlogmas_b	iniSDlogmas_c	iniSDvul_b	iniSDvul_c	iniSDgen_c	");
  fprintf(data1,"finSDlogmas_b	finSDlogmas_c	finSDvul_b	finSDvul_c	finSDgen_c	");
  fprintf(data1,"iniL	finL	initop_c	fintop_c	initop_b	fintop_b	iniOmn	finOmn	");
  fprintf(data1,"alg_mean_m_p	alg_mean_m_c	geom_mean_m_p	geom_mean_m_c	");
  fprintf(data1,"interc_b	expon_b	R2_b	interc_c	expon_c	R2_c	interc	expon	R2	");
  fprintf(data1,"cv_flag	f_herbiv	f_pred	\n"); 

  fclose(data1);

  return;
}

//  ********** invoke the ode solver, i.e., simulate the dynamics of one given system (and calculate mean and SD of biomasses) **********
static void solve_ode(double B[], double meanB[], double CV[], double params[])
{
  int i,flag,ERR_FLAG;
  double t = 0;                            // starting time
  double tout = 0.0001;                        // initial step size
  double *meanB2=(double *) calloc(S+N,sizeof(double));        // mean biomass densities (for comparison)

  gsl_vector_view Bvec = gsl_vector_view_array(B,S+N);        // biomass densities (including nutrients)
  gsl_vector_view Mvec = gsl_vector_view_array(meanB,S+N);      // mean biomass denisties (including resources)
  gsl_vector_view Mvec2 = gsl_vector_view_array(meanB2,S+N);      // mean biomass denisties (including resources), for comparison

  FILE *timeseries;
  if(TIMESERIES)
    Prepare_timeseries_file(timeseries);

  N_Vector y = N_VNew_Serial(S+N);
  N_VSetArrayPointer(B, y);  
  
  void *cvode_mem = CVodeCreate(CV_BDF, CV_NEWTON); 
  CVodeSetUserData(cvode_mem,params);
  CVodeInit(cvode_mem, Cvode_rates, t, y);
  CVodeSStolerances(cvode_mem, eps_rel, eps_abs);
  CVDense(cvode_mem, S+N);
  
// ********* integrate ODE to reach an attractor ********
  ERR_FLAG = 0;
  while (t < tend) {
    flag = CVode(cvode_mem, tout, y, &t, CV_ONE_STEP);

    if (flag != CV_SUCCESS){
    ERR_FLAG = 1;
      break;
    }
      
    if (TIMESERIES)
      Write_timeseries_to_file(B,t,timeseries);
  }

  if(ERR_FLAG == 1){
    for(i=0; i<S; i++){
      if(B[i] < EXTINCT)        // if an error occured, set small biomass densities to 0 
        B[i] = 0;
    }
    t = 0;
    
    while (t < tend) {
      flag = CVode(cvode_mem, tout, y, &t, CV_ONE_STEP);

      if (flag != CV_SUCCESS){
      printf("still a problem...\n");
        break;
      }
      
      if (TIMESERIES)
        Write_timeseries_to_file(B,t,timeseries);
    }
  }

// ********* integrate ODE to calculate mean biomass densities and coefficients of variation ********
  t = tend;
  tout = t + Delta_t;
  while (t < tend + (1+VAR_COEFF)*teval) {
    flag = CVode(cvode_mem, tout, y, &t, CV_NORMAL);

    if(flag != CV_SUCCESS)
      break;
    else
      tout += Delta_t;
   
    if (t < tend + teval)
      gsl_vector_add(&Mvec.vector,&Bvec.vector);

    if (fabs(t - (tend + teval)) < 0.1*Delta_t)
      gsl_vector_scale(&Mvec.vector,Delta_t/teval);

    if (t >= tend + teval) {
      for (i=0; i < S+2; i++)
        CV[i] += (B[i]-meanB[i])*(B[i]-meanB[i]) * Delta_t/teval;
      gsl_vector_add(&Mvec2.vector,&Bvec.vector);          // compute mean biomass densities again
    }
  }    
  
  gsl_vector_scale(&Mvec2.vector,Delta_t/teval);

  double norm1 = gsl_blas_dasum(&Mvec.vector);
  double norm2 = gsl_blas_dasum(&Mvec2.vector);
  
  cv_flag = norm1 - norm2;
  
  Extinction(meanB,S+N);

  for(i=0; i < S+N; i++)
  {
    if(meanB[i] > EXTINCT)
      CV[i] = sqrt(CV[i])/meanB[i];  
    else
      CV[i]=0;
  }

  CVodeFree(&cvode_mem);    
  free(meanB2);

  return;
}

//  ********** set biomass densities below the extinction threshold to 0 **********
static void Extinction(double B[], int N) {
  int i;
  for(i=0; i<N; i++)
    if(B[i] < EXTINCT)
      B[i] = 0;
      
  return;
}

//  ********** this function is required by the cvode-solver **********
static int Cvode_rates(double t, N_Vector y, N_Vector ydot, void *params) {  
  dynamics(N_VGetArrayPointer(y),N_VGetArrayPointer(ydot),params);
  return(0);
}

//  ********** helper function that prints a matrix in the standard output (e.g. a terminal) **********
static void show_matrix(gsl_matrix *A, int N) {
  int i,j;
  for(i = 0; i<N; i++) {
    for(j = 0; j<N; j++)
      printf("%g ",gsl_matrix_get(A,i,j));
    printf("\n");
  }
  return;
}

//  ********** generate output **********
static void output(gsl_matrix *Ap, gsl_vector *mas, double B[], double iniB[], double meanB[], double VarCoeffB[], void *params, FILE *data1) {
  int i,j,k;
  gsl_vector *defvec=gsl_vector_calloc(S);   // default vector: use it, but don't store important data in it!
  gsl_matrix *defmat=gsl_matrix_calloc(S,S); // default matrix: same as with defvec....
  double *pparams=(double *)params;          // pointer to first element of params...

//  ********** persistence, mean biomass densities, and coefficients of variation **********
  double persistence=0,persistence_b=0,persistence_c=0,persistence_n=0;
  double meanB_n=0,meanB_b=0,meanB_c=0,meanB_total=0;
  double CV_n=0,CV_b=0,CV_c=0,CV=0;  

  for(i=0; i<S+N; i++) {
    if(B[i] >= EXTINCT)
    {
      if(i < S) {                            // all species (excluding nutrients)
        persistence++;
        meanB_total+=meanB[i];
        CV+=VarCoeffB[i];
      }
      if(i<S_b) {                            // basal species
        persistence_b++;
        meanB_b+=meanB[i];
        CV_b+=VarCoeffB[i];
      }
      if(i >= S_b && i < S) {                // consumers
        persistence_c++;
        meanB_c+=meanB[i];
        CV_c+=VarCoeffB[i];
      }
      if(i >= S) {                           // nutrients
        persistence_n++;
        meanB_n+=meanB[i];
        CV_n+=VarCoeffB[i];
      }
    }
  }

  if(persistence!=0) {
    meanB_total=meanB_total/persistence;
    CV=CV/persistence;
  }
  else {
    meanB_total=0;
    CV=0;
  }
  if(persistence_b!=0) {
    meanB_b=meanB_b/persistence_b;
    CV_b=CV_b/persistence_b;
  }
  else {
    meanB_b=0;
    CV_b=0;
  }
  if(persistence_c!=0) {
    meanB_c=meanB_c/persistence_c;
    CV_c=CV_c/persistence_c;
  }
  else {
    meanB_c=0;
    CV_c=0;
  }
  if(persistence_n!=0) {
    meanB_n=meanB_n/persistence_n;
    CV_n=CV_n/persistence_n;
  }
  else {
    meanB_n=0;
    CV_n=0;
  }

  persistence=persistence/S;
  persistence_b=persistence_b/S_b;
  persistence_c=persistence_c/S_c;
  persistence_n=persistence_n/N;

//  ********** for final network properties: remove extinct species from the network **********
  gsl_matrix_memcpy(defmat,Ap);
  gsl_vector_memcpy(defvec,mas);
  
  for(i=0;i<S;i++) {
    if(B[i]<1e-6) {
      gsl_vector_set(defvec,i,0);
      for(j=0;j<S;j++) {
        gsl_matrix_set(defmat,i,j,0);
        gsl_matrix_set(defmat,j,i,0);
      }
    }
  }

//  ********** calculate trophic levels **********
  double *TLmin=(double *) calloc(S,sizeof(double));      // memory for minimal trophic levels
  double *TLav=(double *) calloc(S,sizeof(double));       // memory for average trophic levels
  double *TLeff=(double *) calloc(S,sizeof(double));      // memory for effective trophic levels

  double meanTLmin_ini=0,meanTLmin_fin=0;
  double meanTLav_ini=0,meanTLav_fin=0;
  double meanTLeff=0;
  int iniOmn=0, finOmn=0;
  
  Min_TL(Ap,TLmin);                                       // initiale minimale trophische Level berechnen
  Average_TL(Ap,TLav,TLmin);                              // initiale mittlere trophische Level berechnen

  for(i=0;i<S_c;i++) {
    meanTLmin_ini+=TLmin[S_b+i];
    meanTLav_ini+=TLav[S_b+i];
  }
  meanTLmin_ini=meanTLmin_ini/S_c;                        // mean initial minimal trophic level of consumers
  meanTLav_ini=meanTLav_ini/S_c;                          // mean initial average trophic level of consumers


  for(i = 0; i < S_c; i++) {
    k=0;
    for(j = 0; j < S; j++) {    
      if(gsl_matrix_get(Ap,S_b+i,j) == 1 && TLmin[j] >= TLmin[S_b+i])
        k=1;                    
    }
    iniOmn+=k;                                            // initial number of omnivores
  }

  Min_TL(defmat,TLmin);                                   // calculate final minimal trophic levels
  Average_TL(defmat,TLav,TLmin);                          // calculate final average trophic levels
  Eff_TL(defmat,meanB,params,TLeff,TLmin);                // calculate effective trophic levels

  for(i=0;i<S_c;i++) {
    if(B[S_b+i]>=1e-6) {
      meanTLmin_fin+=TLmin[S_b+i];
      meanTLav_fin+=TLav[S_b+i];
      meanTLeff+=TLeff[S_b+i];
    }
  }
  if(persistence_c!=0) {
    meanTLmin_fin=meanTLmin_fin/(S_c*persistence_c);       // mean final minimal trophic level of consumers
    meanTLav_fin=meanTLav_fin/(S_c*persistence_c);         // mean final average trophic level of consumers
    meanTLeff=meanTLeff/(S_c*persistence_c);               // mean effective trophic level of consumers
  }
  else {
    meanTLmin_fin=0;
    meanTLav_fin=0;
    meanTLeff=0;
  }
  
  for(i = 0; i < S_c; i++) {
    k = 0;
    for(j = 0; j < S; j++) {    
      if(gsl_matrix_get(defmat,S_b+i,j) == 1 && TLmin[j] >= TLmin[S_b+i]) 
        k = 1;
    }
    finOmn += k;                                           // final number of omnivores  
  }
    
//  links, body masses, generality, and vulnerability
  double *logmas=(double *) calloc(S,sizeof(double));
  double *vul=(double *) calloc(S,sizeof(double));  
  double *gen=(double *) calloc(S,sizeof(double));
  
  double inivul_b=0,inivul_c=0,inigen_c=0;
  double finvul_b=0,finvul_c=0,fingen_c=0;
  double iniSDvul_b=0,iniSDvul_c=0,iniSDgen_c=0;
  double finSDvul_b=0,finSDvul_c=0,finSDgen_c=0;
  double inilogmas_b=0,inilogmas_c=0,finlogmas_b=0,finlogmas_c=0;
  double iniSDlogmas_b=0,iniSDlogmas_c=0,finSDlogmas_b=0,finSDlogmas_c=0;
  double inidens_c=0,inidens_b=0,findens_c=0,findens_b=0;
  int ini_l=0, fin_l=0, initop_c=0, fintop_c=0, initop_b=0, fintop_b=0;
  
  for(i=0; i<S; i++) {
    for(j=0; j<S; j++) {
      ini_l+=gsl_matrix_get(Ap,i,j);                       // initial number of links
      fin_l+=gsl_matrix_get(defmat,i,j);                   // final number of links
    }
  }

  for(i=0;i<S;i++) {
    logmas[i]=log10(gsl_vector_get(mas,i));
  
    gsl_vector_view tempp_in=gsl_matrix_row(Ap,i);
    gen[i]=gsl_blas_dasum(&tempp_in.vector);               // generality

    gsl_vector_view tempp_out=gsl_matrix_column(Ap,i);
    vul[i]=gsl_blas_dasum(&tempp_out.vector);              // vulnerability (out-degree distributions)
  }
  
  for(i=0;i<S_b;i++) {
    inivul_b+=vul[i]/S_b;                                  // mean initial vulnerability of basals
    inilogmas_b+=logmas[i]/S_b;                            // mean initial log body mass of basals
    inidens_b+= iniB[i]/gsl_vector_get(mas,i);             // initial total population density of basals 
  }
  for(i=0;i<S_c;i++) {
    inivul_c+=vul[S_b+i]/S_c;                              // mean initial vulnerability of consumers
    inigen_c+=gen[S_b+i]/S_c;                              // mean initial generality of consumers
    inilogmas_c+=logmas[S_b+i]/S_c;                        // mean initial log body mass of consumers
    inidens_c+= iniB[S_b+i]/gsl_vector_get(mas,S_b+i);     // initial total populations density of all consumers
  }
  
  for(i=0;i<S_c;i++) {
      if(vul[S_b+i] == 0 && gen[S_b+i] > 0) 
      initop_c+=1;                                         // initial number of top predators
  }
  
  for(i=0;i<S_b;i++) {
      if(vul[i] == 0) 
      initop_b+=1;                                         // initial number of predator-free basals 
  }
  
  for(i=0;i<S_b;i++) {
    iniSDvul_b+=(vul[i]-inivul_b)*(vul[i]-inivul_b)/S_b;
    iniSDlogmas_b+=(logmas[i]-inilogmas_b)*(logmas[i]-inilogmas_b)/S_b;
  }
  for(i=0;i<S_c;i++) {
    iniSDvul_c+=(vul[S_b+i]-inivul_c)*(vul[S_b+i]-inivul_c)/S_c;
    iniSDgen_c+=(gen[S_b+i]-inigen_c)*(gen[S_b+i]-inigen_c)/S_c;
    iniSDlogmas_c+=(logmas[S_b+i]-inilogmas_c)*(logmas[S_b+i]-inilogmas_c)/S_c;
  }
  iniSDvul_b=sqrt(iniSDvul_b);                             // SD initial vulnerability of basals
  iniSDlogmas_b=sqrt(iniSDlogmas_b);                       // SD initial log body mass of basals
  iniSDvul_c=sqrt(iniSDvul_c);                             // SD initial vulnerability of consumers
  iniSDgen_c=sqrt(iniSDgen_c);                             // SD initial generality of consumers
  iniSDlogmas_c=sqrt(iniSDlogmas_c);                       // SD initial log body mass of consumers

//  final quantities:
  for(i=0;i<S;i++) {
    gsl_vector_view tempp_in=gsl_matrix_row(defmat,i);
    gen[i]=gsl_blas_dasum(&tempp_in.vector);                // generality

    gsl_vector_view tempp_out=gsl_matrix_column(defmat,i);
    vul[i]=gsl_blas_dasum(&tempp_out.vector);               // vulnerability (out-degree distribution)
  }

  for(i=0;i<S_b;i++) {
    if(B[i]>=1e-6) {
      finvul_b+=vul[i]/(S_b*persistence_b);                 // mean final vulnerability of basals
      finlogmas_b+=logmas[i]/(S_b*persistence_b);           // mean final log body mass of basals
      findens_b+= meanB[i]/gsl_vector_get(mas,i);           // final total population density of basals
    }
  }
  for(i=0;i<S_c;i++) {
    if(B[S_b+i]>=1e-6) {
      finvul_c+=vul[S_b+i]/(S_c*persistence_c);             // mean fianl vulnerability of consumers
      fingen_c+=gen[S_b+i]/(S_c*persistence_c);             // mean final genarality of consumers
      finlogmas_c+=logmas[S_b+i]/(S_c*persistence_c);       // mean final log body mass of consumers
      findens_c+= meanB[S_b+i]/gsl_vector_get(mas,S_b+i);   // final total population density of all consumers
    }
  }
  
  for(i=0;i<S_c;i++) {
    if(B[i]>=1e-6)
      if(vul[S_b+i] == 0 && gen[S_b+i] > 0) 
        fintop_c+=1;                                        // final number of top predators
  }
  
  for(i=0;i<S_b;i++) {
    if(B[i]>=1e-6)
      if(vul[i] == 0) 
        fintop_b+=1;                                        // final number of predator-free basal species 
  }
  
  for(i=0;i<S_b;i++) {
    if(B[i]>=1e-6) {
      finSDvul_b+=(vul[i]-finvul_b)*(vul[i]-finvul_b)/(S_b*persistence_b);
      finSDlogmas_b+=(logmas[i]-finlogmas_b)*(logmas[i]-finlogmas_b)/(S_b*persistence_b);
    }
  }
  for(i=0;i<S_c;i++) {
    if(B[S_b+i]>=1e-6) {
      finSDvul_c+=(vul[S_b+i]-finvul_c)*(vul[S_b+i]-finvul_c)/(S_c*persistence_c);
      finSDgen_c+=(gen[S_b+i]-fingen_c)*(gen[S_b+i]-fingen_c)/(S_c*persistence_c);
      finSDlogmas_c+=(logmas[S_b+i]-finlogmas_c)*(logmas[S_b+i]-finlogmas_c)/(S_c*persistence_c);
    }
  }
  finSDvul_b=sqrt(finSDvul_b);                              // SD final vulnerability of basals
  finSDlogmas_b=sqrt(finSDlogmas_b);                        // SD final log body mass of basals
  finSDvul_c=sqrt(finSDvul_c);                              // SD final vulnerability of consumers
  finSDgen_c=sqrt(finSDgen_c);                              // SD final generality of consumers
  finSDlogmas_c=sqrt(finSDlogmas_c);                        // SD final log body mass of consumers


//  functional diversity
  double initial_funct_div,final_funct_div;                 // initial and final functional diversity
  initial_funct_div=functional_diversity(mas);
  final_funct_div=functional_diversity(defvec);

//  biomass flows
  double *energyflows=(double *) calloc(6,sizeof(double));  // 0: nutrient influx, 1: nutrient-basal, 2: basal-consumer, 3: igp, 4: basal resp., 5: cons. resp.
  flows(meanB,params,energyflows);

//  mean individual body masses (alg. mean basal, alg. mean consumer, geom. mean basal, geom. mean consumer)
  double *mean_body_masses = (double *) calloc(4,sizeof(double));  
  mean_body_mass(meanB, mas, mean_body_masses);

//  mass-abundance relation
  double *fitresults=(double *) calloc(9,sizeof(double));
  mass_abundance(meanB,mas,fitresults);

//  write data to file
  fprintf(data1,"%d	%d	%d	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	",Web_Id, S_b, S_c, hill, C_interference, a_c, a_p, h_c, h_p, persistence, persistence_b, persistence_c, persistence_n);
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	",meanB_total, meanB_b, meanB_c, meanB_n, CV, CV_b, CV_c, CV_n);
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	",initial_funct_div, final_funct_div, inidens_b, inidens_c, findens_b, findens_c);
  for(i=0;i<N;i++)
    fprintf(data1,"%.6g	",*(pparams+2*S*S+5*S+2*S_b*N+i));	// nutrient supply concentrations
  for(i=0;i<6;i++)
    fprintf(data1,"%.6g	",energyflows[i]);
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	",meanTLmin_ini,meanTLmin_fin,meanTLav_ini,meanTLav_fin,meanTLeff);
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	",inilogmas_b,inilogmas_c,inivul_b,inivul_c,inigen_c);	
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	",finlogmas_b,finlogmas_c,finvul_b,finvul_c,fingen_c);
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	",iniSDlogmas_b,iniSDlogmas_c,iniSDvul_b,iniSDvul_c,iniSDgen_c);
  fprintf(data1,"%.6g	%.6g	%.6g	%.6g	%.6g	",finSDlogmas_b,finSDlogmas_c,finSDvul_b,finSDvul_c,finSDgen_c);
  fprintf(data1,"%d	%d	%d	%d	%d	%d	%d	%d	", ini_l, fin_l, initop_c, fintop_c, initop_b, fintop_b, iniOmn, finOmn);
  for(i=0; i<4; i++)
    fprintf(data1,"%.6g	",mean_body_masses[i]);
  for(i=0;i<9;i++)
    fprintf(data1,"%.6g	",fitresults[i]);
  fprintf(data1,"%.6g	%.6g	%.6g\n", cv_flag,f_herbiv,f_pred);

  free(TLmin);
  free(TLav);
  free(TLeff);
  free(vul);
  free(gen);
  free(energyflows);
  free(fitresults);
  free(mean_body_masses);
  gsl_vector_free(defvec);
  gsl_matrix_free(defmat);

  return;
}

//  ********** calculate algebraic and geometric mean body masses **********
static void mean_body_mass(double meanB[], gsl_vector *mass, double mean_body_masses[]) {
  int i;
  double temp1,temp2,temp3,temp4;
  double *populations = (double *) calloc(S,sizeof(double));
  
  for(i = 0; i<S; i++)
   *(populations + i) = meanB[i] / gsl_vector_get(mass,i);
  
  temp1 = 0;
  temp2 = 0;
  for(i = 0; i< S_b; i++) {
    temp1 += meanB[i];
    temp2 += populations[i];
  }
  mean_body_masses[0] = temp1/temp2;                          // algebraic mean body mass of basal individuals

  temp1 = 0;
  temp2 = 0;
  for(i = S_b; i < S; i++) {
    temp1 += meanB[i];
    temp2 += populations[i];
  }
  mean_body_masses[1] = temp1/temp2;                          // algebraic mean body mass of consumer individuals
  
  temp1 = 0;
  temp2 = 0;
  for(i = 0; i< S_b; i++) {
    temp1 += populations[i] * log(gsl_vector_get(mass,i));
    temp2 += populations[i];
  }
  mean_body_masses[2] = exp(temp1 / temp2);                   // geometric mean body mass of basal individuals
  
  temp1 = 0;
  temp2 = 0;
  for(i = S_b; i< S; i++) {
    temp1 += populations[i] * log(gsl_vector_get(mass,i));
    temp2 += populations[i];
  }
  mean_body_masses[3] = exp(temp1 / temp2);                   // geometric mean body mass of consumer individuals
  
  free(populations);
  
  return;
}  

//  ********** calculate energy flows **********
static void flows(double meanB[], void *params, double energyflows[]) {
  int i,j;

  gsl_vector_view TB_vec = gsl_vector_view_array(meanB,S+N);  // vector with all biomass densities and nutrient concentration
  gsl_vector *TBvec = &TB_vec.vector;

  gsl_vector *Ivec=gsl_vector_calloc(S);                      // ingestion: nutrient uptake, consumption, predation 
  gsl_vector *Dvec=gsl_vector_calloc(S);                      // death (predation losses)
  gsl_vector *Xvec=gsl_vector_calloc(S);                      // respiration
  gsl_vector *N_in=gsl_vector_calloc(N);                      // influx of nutrients
  gsl_vector *N_out=gsl_vector_calloc(N);                     // outflux of nutrients
  gsl_vector *N_up=gsl_vector_calloc(N);                      // uptake of nutrients by plants

  net_rates(TBvec,Ivec,Dvec,Xvec,N_in,N_out,N_up,params);     // calculate net rates

  energyflows[0] = gsl_blas_dasum(N_in);                      // flux 0: nutrient influx 
  energyflows[1] = gsl_blas_dasum(N_up);                      // flux 1: nutrient uptake by basal species (plants) 

  for(i=0;i<S_b;i++)
    energyflows[2] += gsl_vector_get(Dvec,i);                 // flux 2: basal to consumers

  for(i=0;i<S_c;i++)
    energyflows[3] += gsl_vector_get(Dvec,S_b+i);             // flux 3: consumers intraguild

  for(i=0;i<S_b;i++)
    energyflows[4] += gsl_vector_get(Xvec,i);                 // flux 4: basal respiration

  for(i=0;i<S_c;i++)
    energyflows[5] += gsl_vector_get(Xvec,S_b+i);             // flux 5: consumer respiration

//  ********** free memory ************
  gsl_vector_free(Ivec);
  gsl_vector_free(Dvec);
  gsl_vector_free(Xvec);
  gsl_vector_free(N_in);
  gsl_vector_free(N_out);
  gsl_vector_free(N_up);

  return;
}

//  ********** determine mass-abundance relation **********
static void mass_abundance(double meanB[], gsl_vector *mas, double fitresults[]) {
  int i,j;
  int c=0,b=0;

  double exponent,exponent_b,exponent_c,intercept,intercept_b,intercept_c,R2,R2_b,R2_c;

  for(i=0;i<S_b;i++) {
    if(meanB[i]>0)
      b++;
  }
  for(i=0;i<S_c;i++) {
    if(meanB[i+S_b]>0)
      c++;
  }

  double *logmas_b=(double *) calloc(b,sizeof(double));
  double *logmas_c=(double *) calloc(c,sizeof(double));
  double *logmas=(double *) calloc(b+c,sizeof(double));
  double *logabundance_b=(double *) calloc(b,sizeof(double));
  double *logabundance_c=(double *) calloc(c,sizeof(double));
  double *logabundance=(double *) calloc(b+c,sizeof(double));

//  data: log body masses and biomass densities of basal species and consumers
  j=0;
  for(i=0;i<S_b;i++) {
    if(meanB[i]>0) {
      logmas_b[j]=log(gsl_vector_get(mas,i));        // basal species
      logabundance_b[j]=log(meanB[i]);
      j++;
    }
  }
  j=0;
  for(i=0;i<S_c;i++) {
    if(meanB[i+S_b]>0) {
      logmas_c[j]=log(gsl_vector_get(mas,i+S_b));    // consumers
      logabundance_c[j]=log(meanB[i+S_b]);
      j++;
    }
  }
  for(i=0;i<b;i++) {
    logmas[i]=logmas_b[i];
    logabundance[i]=logabundance_b[i];               // basal species and consumers together
  }
  for(i=0;i<c;i++) {
    logmas[b+i]=logmas_c[i];
    logabundance[b+i]=logabundance_c[i];
  }

//  fit a linear model to the data
  double meanX=0,meanY=0,SSXY=0,SSX=0,SS_tot=0,SS_err=0;
  for(i=0;i<b;i++) {                                 // first the basal species...
    meanX+=logmas_b[i]/b;
    meanY+=logabundance_b[i]/b;
  }
  for(i=0;i<b;i++) {
    SSXY+=(logmas_b[i]-meanX)*(logabundance_b[i]-meanY);
    SSX+=(logmas_b[i]-meanX)*(logmas_b[i]-meanX);
  }
  exponent_b=SSXY/SSX;
  intercept_b=meanY-exponent_b*meanX;

  for(i=0;i<b;i++) {                                 // calculate r^2
    SS_tot+=(logabundance_b[i]-meanY)*(logabundance_b[i]-meanY);
    SS_err+=(logabundance_b[i]-intercept_b-exponent_b*logmas_b[i])*(logabundance_b[i]-intercept_b-exponent_b*logmas_b[i]);
  }
  R2_b=1-SS_err/SS_tot;  

  meanX=0;
  meanY=0;
  SSXY=0;
  SSX=0;
  SS_tot=0;
  SS_err=0;

  for(i=0;i<c;i++)                                   // now the consumers...
  {
    meanX+=logmas_c[i]/c;
    meanY+=logabundance_c[i]/c;
  }
  for(i=0;i<c;i++) {
    SSXY+=(logmas_c[i]-meanX)*(logabundance_c[i]-meanY);
    SSX+=(logmas_c[i]-meanX)*(logmas_c[i]-meanX);
  }
  exponent_c=SSXY/SSX;
  intercept_c=meanY-exponent_c*meanX;

  for(i=0;i<c;i++) {
    SS_tot+=(logabundance_c[i]-meanY)*(logabundance_c[i]-meanY);
    SS_err+=(logabundance_c[i]-intercept_c-exponent_c*logmas_c[i])*(logabundance_c[i]-intercept_c-exponent_c*logmas_c[i]);
  }
  R2_c=1-SS_err/SS_tot;  

  meanX=0;
  meanY=0;
  SSXY=0;
  SSX=0;
  SS_tot=0;
  SS_err=0;

  for(i=0;i<b+c;i++) {                               // ... and finally all species together
    meanX+=logmas[i]/(b+c);
    meanY+=logabundance[i]/(b+c);
  }
  for(i=0;i<b+c;i++) {
    SSXY+=(logmas[i]-meanX)*(logabundance[i]-meanY);
    SSX+=(logmas[i]-meanX)*(logmas[i]-meanX);
  }
  exponent=SSXY/SSX;
  intercept=meanY-exponent*meanX;

  for(i=0;i<b+c;i++) {
    SS_tot+=(logabundance[i]-meanY)*(logabundance[i]-meanY);
    SS_err+=(logabundance[i]-intercept-exponent*logmas[i])*(logabundance[i]-intercept-exponent*logmas[i]);
  }
  R2=1-SS_err/SS_tot;  

  fitresults[0]=intercept_b;
  fitresults[1]=exponent_b;
  fitresults[2]=R2_b;
  fitresults[3]=intercept_c;
  fitresults[4]=exponent_c;
  fitresults[5]=R2_c;
  fitresults[6]=intercept;
  fitresults[7]=exponent;
  fitresults[8]=R2;

  free(logmas_b);
  free(logmas_c);
  free(logmas);
  free(logabundance_b);
  free(logabundance_c);
  free(logabundance);
  
  return;
}

//  ********** calculate minimal distance of the species from the nutrients **********
static void Min_TL(gsl_matrix *Ap,double TLmin[]) {
  int i,j,k,counter;

  for(i=0;i<S;i++) {
    gsl_vector_view tempp=gsl_matrix_row(Ap,i);  
    if(gsl_vector_isnull(&tempp.vector)==1)
      TLmin[i]=1;                                  // determine species on the first trophic level
    else
      TLmin[i]=0;
  }  
  k=2;
  while(k<=S) {
    counter=0;
    for(i=0;i<=S-1;i++) {
      if(TLmin[i]==0) {
        for(j=0;j<=S-1;j++) {
          if(TLmin[j]==k-1&&gsl_matrix_get(Ap,i,j)!=0) {  
            TLmin[i]=k;                            // determine species on level k >= 2
            counter=counter+1;
          }
        }
      }
    }
    if(counter==0)
      k=S;  
    k=k+1;
  }  

  return;
}

//  ********** calculate mean distance of the species from the nutrients **********
static void Average_TL(gsl_matrix *Ap, double TLav[], double TLmin[]) {
// idea: solve an equation of the form (matrix)Ta * (vector)xa = (vector)b 
  int i,j;
  double norm;
  int sig=1;
  gsl_vector *b=gsl_vector_calloc(S);            // inhomogeneity 
  gsl_vector *xa=gsl_vector_calloc(S);           // unknown variables
  gsl_permutation *p=gsl_permutation_calloc(S);  // permutation vector
  gsl_matrix *Ta=gsl_matrix_calloc(S,S);         // matrix of the system of equations

  gsl_matrix_memcpy(Ta,Ap);

  for(i=0;i<S;i++) {
    gsl_vector_view tempp=gsl_matrix_row(Ta,i);
    norm=gsl_blas_dasum(&tempp.vector);
    
    if(norm!=0) {
      for(j=0;j<S;j++)
        gsl_matrix_set(Ta,i,j,gsl_matrix_get(Ta,i,j)/norm);
    }

    if(TLmin[i]==0) {
      for(j=0;j<S;j++)  
        gsl_matrix_set(Ta,i,j,0);
    }
  }
  for(i=0;i<S;i++)
    gsl_matrix_set(Ta,i,i,gsl_matrix_get(Ta,i,i)-1);

  for(i=0;i<S;i++) {
    if(TLmin[i]!=0)
      gsl_vector_set(b,i,-1);
  }

  gsl_linalg_LU_decomp(Ta,p,&sig);               // decompose matrix Ta... 
  gsl_linalg_LU_solve(Ta,p,b,xa);                // ... and solve the system

  for(i=0;i<S;i++)
    TLav[i]=gsl_vector_get(xa,i);

  gsl_matrix_free(Ta);
  gsl_vector_free(xa);
  gsl_permutation_free(p);
  gsl_vector_free(b);

  return;
}

//  ********** calculate effective distance of the species from the nutrients **********
static void Eff_TL(gsl_matrix *Ap, double meanB[], void *params, double TLeff[], double TLmin[]) {
  int i,j;
  double *pparams=(double *)params;                                   // pointer to the first element of params

  gsl_vector *Bvec=gsl_vector_calloc(S);                              // mean biomasses  
  gsl_vector_const_view B_vec1=gsl_vector_const_view_array(meanB,S);
  gsl_vector_memcpy(Bvec,&B_vec1.vector);

  gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);           // a_ij*f_ij
  gsl_matrix *Amat=&A_mat.matrix;

  gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);       // H_ij
  gsl_matrix *Hmat=&H_mat.matrix;

  gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);   // interference vector
  gsl_vector *Cvec=&C_vec.vector;

  gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);   // hill coefficients
  gsl_vector *Hvec=&H_vec.vector;  

  gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);   // body masses
  gsl_vector *Mvec=&M_vec.vector;  

  gsl_vector *Prey=gsl_vector_calloc(S);                              // prey biomasses

  gsl_matrix *Effmat=gsl_matrix_calloc(S,S);
  gsl_matrix *mmat=gsl_matrix_calloc(S,S);
  gsl_vector *rvec=gsl_vector_calloc(S);
  gsl_vector *svec=gsl_vector_calloc(S);

//  ************ in case of Type 3 functional response ***********
  gsl_vector_memcpy(Prey,Bvec);

  for(i=0;i<S;i++) {
    if(gsl_vector_get(Prey,i)<0)
      gsl_vector_set(Prey,i,0);                       // paranoia (this doesn't mean that you can safely remove these lines!)

    gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));
  }

//  ************* functional responses zusammenbauen *************
  gsl_matrix_memcpy(mmat,Hmat);
  gsl_matrix_mul_elements(mmat,Amat);

  gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);    // handling time term: r_i=sum_j a_ij*h_ij*Prey_j
  gsl_vector_memcpy(svec,Bvec);
  gsl_vector_mul(svec,Cvec);                          // predator interference: s_i=C_i*B_i
  gsl_vector_add(rvec,svec);
  gsl_vector_add_constant(rvec,1);
  gsl_vector_mul(rvec,Mvec);                          // rvec: denominator of functional response

  for(i=0;i<S;i++) {
    for(j=0;j<S;j++)
      gsl_matrix_set(Effmat,i,j,gsl_matrix_get(Amat,i,j)*gsl_vector_get(Prey,j)/gsl_vector_get(rvec,i));
  }

  Average_TL(Effmat, TLeff, TLmin);

//  ********** free memory ************
  gsl_vector_free(Bvec);
  gsl_vector_free(Prey);
  gsl_vector_free(svec);
  gsl_vector_free(rvec);
  gsl_matrix_free(mmat);
  gsl_matrix_free(Effmat);

  return;
}

//  ********** calculate functional diversity **********
double functional_diversity(gsl_vector *mas) {
  int i,j;
  double max;                    // maximum of attack rates
  double M;                      // body mass of a species
  double M_prey;                 // body mass of a prey species
  double attackrate;             // atack rate of a species 

  double L;                      // log_10 length of body-mass axis

  if(UNIFORM_MASSES == 0)
    L=10+fabs(zeta_c - zeta_b);
  else
    L=10 + fabs(m_a_min - m_p_min);

  double N=1200;                 // number of supporting points where the integral is evaluated
  double l=L/N;                  // length of a single interval
  double s;                      // point at which the function is evaluated
  double I=0;                    // value of the integral (measure of functional diversity)

  for(j=0;j<N;j++) {    
    max=0;
    s=(j+0.5)*l;
    M_prey=pow(10,s);

    for(i=0;i<S_c;i++) {
      M=gsl_vector_get(mas,i+S_b);
      attackrate=pow(M/(M_prey*R_opt)*exp(1-(M/(M_prey*R_opt))),g);

      if(attackrate<cutoff)
        attackrate=0;

      if(attackrate>max)
        max=attackrate;
    }

    I+=max*l;
  }

  double correction = 2;
  
  for(j=0; j<1000; j++) {        // calculate correction term
    s = (j+0.5) * 0.002;
    M_prey = pow(10,s);
    
    attackrate = pow(100./(M_prey*R_opt)*exp(1.-(100./(M_prey*R_opt))),g);

    if(attackrate<cutoff)
      attackrate=0;
      
    correction -= attackrate*0.002;
  }
  
  I = I/(L-correction);          // normalisation

  return(I);
}

//  ********** define population dynamics equations **********
static void dynamics(double B[], double Bdot[], void *params) {
  int i;

  gsl_vector_view TB_vec = gsl_vector_view_array(B,S+N);      // vector with all biomass densities and nutrient concentration
  gsl_vector *TBvec = &TB_vec.vector;

  gsl_vector_view Bdot_vec=gsl_vector_view_array(Bdot,S);     // vector for time derivatives of biomass densities
  gsl_vector *Bdotvec=&Bdot_vec.vector;

  gsl_vector_view Ndot_vec=gsl_vector_view_array(Bdot+S,N);   // vector for time derivatives of nutrient concentrations
  gsl_vector *Ndotvec=&Ndot_vec.vector;

  gsl_vector *Ivec=gsl_vector_calloc(S);                      // ingestion: nutrient uptake, consumption, predation 
  gsl_vector *Dvec=gsl_vector_calloc(S);                      // death (predation losses)
  gsl_vector *Xvec=gsl_vector_calloc(S);                      // respiration
  gsl_vector *N_in=gsl_vector_calloc(N);                      // influx of nutrients
  gsl_vector *N_out=gsl_vector_calloc(N);                     // outflux of nutrients
  gsl_vector *N_up=gsl_vector_calloc(N);                      // uptake of nutrients by plants

  net_rates(TBvec,Ivec,Dvec,Xvec,N_in,N_out,N_up,params);     // calculate net rates

//  ************* assmble dynamical equations (Bdot) *************
  gsl_vector_memcpy(Ndotvec,N_in);
  gsl_vector_sub(Ndotvec,N_out);
  gsl_vector_sub(Ndotvec,N_up);

  gsl_vector_memcpy(Bdotvec,Ivec);
  gsl_vector_sub(Bdotvec,Dvec);
  gsl_vector_sub(Bdotvec,Xvec);

//  ************* control of negative biomass densities *************
  for(i=0; i<S; i++) {
    if(B[i] < EXTINCT)
      Bdot[i] = -B[i];
  }

//  ********** free memory ************
  gsl_vector_free(Ivec);
  gsl_vector_free(Dvec);
  gsl_vector_free(Xvec);
  gsl_vector_free(N_in);
  gsl_vector_free(N_out);
  gsl_vector_free(N_up);

  return;
}

//  ********** calculate biomass flow rates between species and resources **********
static void net_rates(gsl_vector *TBvec, gsl_vector *Ivec, gsl_vector *Dvec, gsl_vector *Xvec, 
                gsl_vector *N_in, gsl_vector *N_out, gsl_vector *N_up, void *params) {
  int i;
  double *pparams=(double *)params;                                     // pointer to the first element of params
  
  gsl_vector_view B_vec=gsl_vector_subvector(TBvec,0,S);                // all biomasses (plants + animals)
  gsl_vector *Bvec=&B_vec.vector;

  gsl_vector_view PB_vec=gsl_vector_subvector(TBvec,0,S_b);             // basal biomasses
  gsl_vector *PBvec=&PB_vec.vector;

  gsl_vector_view NB_vec=gsl_vector_subvector(TBvec,S,N);               // nutrient concentrations
  gsl_vector *NBvec=&NB_vec.vector;

  gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);             // a_ij*f_ij
  gsl_matrix *Amat=&A_mat.matrix;

  gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);         // H_ij
  gsl_matrix *Hmat=&H_mat.matrix;

  gsl_vector_view R_vec=gsl_vector_view_array(pparams+2*S*S,S);         // respiration rates
  gsl_vector *Rvec=&R_vec.vector;

  gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);     // interference vector
  gsl_vector *Cvec=&C_vec.vector;

  gsl_vector_view P_vec=gsl_vector_view_array(pparams+2*S*S+2*S,S);     // basal vector (1 for basal species, 0 otherwise)
  gsl_vector *Pvec=&P_vec.vector;

  gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);     // Hill coefficients
  gsl_vector *Hvec=&H_vec.vector;  

  gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);     // body masses
  gsl_vector *Mvec=&M_vec.vector;  

  gsl_vector_view E_vec=gsl_vector_view_array(pparams+2*S*S+5*S,S);     // assimilation efficiencies
  gsl_vector *Evec=&E_vec.vector;  

  gsl_matrix_view K_mat=gsl_matrix_view_array(pparams+2*S*S+6*S,S_b,N); // nutrient uptake half saturation densities
  gsl_matrix *Kmat=&K_mat.matrix;                                       // (P rows, N columns)

  gsl_vector_view S_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N,N);        // nutrient supply concentrations
  gsl_vector *Svec=&S_vec.vector;  

  gsl_vector_view U_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N,S_b);    // max. nutrient uptake rates
  gsl_vector *Uvec=&U_vec.vector;  

  gsl_vector_view Cc_vec=gsl_vector_view_array(pparams+2*S*S+6*S+S_b*N+N+S_b,N); // nutrient contents in plants
  gsl_vector *Ccvec=&Cc_vec.vector;  

  gsl_vector *Prey=gsl_vector_calloc(S);                                // prey biomasses

  gsl_matrix *mmat=gsl_matrix_calloc(S,S);
  gsl_vector *rvec=gsl_vector_calloc(S);
  gsl_vector *svec=gsl_vector_calloc(S);
  gsl_vector *nvec=gsl_vector_calloc(N);
  gsl_vector *mvec=gsl_vector_calloc(N);
  gsl_vector *pvec=gsl_vector_calloc(S_b);

//  ************ in case of Type 3 functional response ***********
  gsl_vector_memcpy(Prey,Bvec);

  for(i=0;i<S;i++) {
    if(gsl_vector_get(Prey,i)<0)
      gsl_vector_set(Prey,i,0);                            // paranoia (this doesn't mean that you can safely remove these lines!)

    gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));
  }

//  ************* functional responses zusammenbauen *************
  gsl_matrix_memcpy(mmat,Hmat);
  gsl_matrix_mul_elements(mmat,Amat);

  gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);         // handling time term: r_i=sum_j a_ij*h_ij*Prey_j
  gsl_vector_memcpy(svec,Bvec);
  gsl_vector_mul(svec,Cvec);                               // predator interference: s_i=C_i*B_i
  gsl_vector_add(rvec,svec);
  gsl_vector_add_constant(rvec,1);
  gsl_vector_mul(rvec,Mvec);                               // rvec: denominator of functional response

  gsl_vector_memcpy(svec,Evec);                            // assimilation efficiency (prey-specific)
  gsl_vector_mul(svec,Prey);
  gsl_blas_dgemv(CblasNoTrans,1,Amat,svec,0,Ivec);         // numerator of prey intake term: I_i = sum_j e_j*A_ij*Prey_j

  gsl_vector_div(Ivec,rvec);                               // divide by denominator...
  gsl_vector_mul(Ivec,Bvec);                               // multiply with target species' biomass

  gsl_vector_memcpy(Xvec,Rvec);                            // compute the total respiration rate:
  gsl_vector_mul(Xvec,Bvec);                               // per unit biomass respiration times biomass density

  gsl_vector_memcpy(svec,Bvec);                            // predator biomass
  gsl_vector_div(svec,rvec);                               // ... divide by denominator of the functional response
  gsl_blas_dgemv(CblasTrans,1,Amat,svec,0,Dvec);           // per capita predation losses:
  gsl_vector_mul(Dvec,Prey);                               // multiply with target species' biomass

//  ************* terms for basal species *************
  gsl_vector_view u_vec=gsl_vector_subvector(rvec,0,S_b);  // basal species growth rates
  gsl_vector *uvec=&u_vec.vector;
  for(i=0;i<S_b;i++) {
    gsl_vector_view K_vec=gsl_matrix_row(Kmat,i);
    gsl_vector *Kvec=&K_vec.vector;                        // Kvec is a vector of length N (Kvec_j=Kmat_ij)

    gsl_vector_memcpy(nvec,NBvec);
    gsl_vector_add(nvec,Kvec);
    gsl_vector_memcpy(mvec,NBvec);
    gsl_vector_div(mvec,nvec);                             // m_j=NB_j/(Kvec_j+NB_j)

    gsl_vector_set(uvec,i,gsl_vector_min(mvec));           // remember: uvec is just a vector_view of rvec!
  }
  gsl_vector_mul(uvec,Uvec);
  gsl_vector_mul(uvec,PBvec);                              // uvec_i=PB_i*mas_i^-0.25*Min_j(NB_j/(NB_j+Kmat_ij))

  gsl_vector_mul(rvec,Pvec);                               // growth function only for plant species
  gsl_vector_add(Ivec,rvec);

//  ********** put the rates for the nutrient dynamics together **********
  gsl_vector_memcpy(N_in,Svec);
  gsl_vector_scale(N_in,D_N);                              // nutrient influx: D_N * S_i
  
  gsl_vector_memcpy(N_out,NBvec);
  gsl_vector_scale(N_out,D_N);                             // nutrient outflux: D_N * N_i
  
  gsl_vector_memcpy(N_up,Ccvec);
  gsl_vector_scale(N_up,gsl_blas_dasum(uvec));             // nutrient uptake: C_i * sum_j (r_j*G_j*P_j)

//  ********** free memory ************
  gsl_vector_free(Prey);
  gsl_vector_free(svec);
  gsl_vector_free(rvec);
  gsl_vector_free(nvec);
  gsl_vector_free(mvec);
  gsl_vector_free(pvec);
  gsl_matrix_free(mmat);
  
  return;
}

//  ********** write the time series of the last run to an output file **********
static void Prepare_timeseries_file(FILE *timeseries) {
  timeseries = fopen(timeseries_file,"w");
  fclose(timeseries);
  
  return;
}

static void Write_timeseries_to_file(double B[], double t, FILE *timeseries) {
  int i;
  timeseries = fopen(timeseries_file,"a");
  
  fprintf(timeseries,"%g ",t);
  for(i=0; i < S+2; i++)
    fprintf(timeseries,"%.8g ", B[i]);
  fprintf(timeseries,"\n");
  
  fclose(timeseries);
  
  return;
}


