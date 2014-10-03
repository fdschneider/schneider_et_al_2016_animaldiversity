/*
            Project: Predator Diversity and Ecosystem Functioning
                                  ------
                           Population Dynamics Code
                              Version 2.0
                      last edit: 06.10.2011 by Christian

        Copyright (C) 2014 Christian Guill & Florian D. Schneider

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
      
    Contact the authors: [Christian Guill](C.P.Guill@uva.nl), [Florian D. Schneider](florian.schneider@univ-montp2.fr)
*/

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>	
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>			// random number generator
#include <gsl/gsl_blas.h>			// linear algebra with vecors and matrices
#include <gsl/gsl_randist.h>			// random number distributions
#include <gsl/gsl_sort_vector.h>		// vector operations
#include <gsl/gsl_errno.h>			// errorhandler
#include <gsl/gsl_linalg.h>			// linear algebra for TL (average)-calculation 
#include <gsl/gsl_permutation.h>		// permutations
#include <gsl/gsl_histogram.h>			// histograms
#include <gsl/gsl_odeiv.h>              // ODE solver 
#include <gsl/gsl_fit.h>			// linear regression

//	basic simulation variables:
int Web_Id=0;
int S;  					// total species number
int S_c=30;					// number of predator species (overwritten)
int S_b=30;					// number of basal species
int N=2;					// number of resources

double tend=150000;  					// time at which simulation runs stop
double teval=10000;						// length of time interval from which to determine mean(B) and SD(B)
int seed=20000;							// seed for random number generator
int Imax=10;							// number of replicates per diversity level
int DivLev=101;							// number of diversity levels
int DivStep=1;							// step size of diversity levels
char filename[256]="pdef_2.0_";			// name of output file

//	Parameter zur Bestimmung der Topologie
double zeta_c=7;  				// log_10 of mean of log-normal distribution of predator body masses 
double sigma_c=4;					// sigma of lognormal distribution (width) of predator body masses 
double cutoff_c=1e5;				// 1/2 range of species body-mass distribution (defining size of smallest and largest potential species) 
double zeta_b=5;					// log_10 of mean of log-normal distribution of basal species body masses 
double sigma_b=4;					// sigma of lognormal distribution (width) of basal species body masses 
double cutoff_b=1e5;				// 1/2 range of species body-mass distribution (defining size of smallest and largest potential species) 

double b1=0.01;				//	Coefficient of opt. bodymass ratio
double b2=1.0;					// exponent of opt. bodymass ratio	
double b3=0.1;					// minimal sigma of link distribution
double b4=0.05;					// exponent of sigma of link distribution
double cutoff2=0.01;			// cutoff for assignment of links

double R_opt=100;			// optimal predator-prey bodymass-ratio
double g=2;					// gamma, breadth of Ricker-function 

//	parameters of functional response
double a_0=200;				// coefficient of attack rate
double a_c=0.25;			// exponent of predator bodymass in attack rate 
double a_p=0.25;			// exponent of prey bodymass in attack rate 
double h_0=0.4;				// coefficient of handling times
double h_c=-0.75;			// exponent of predator bodymass in handling time
double h_p=0;				// exponent of Brey bodymass in handling time
double mu_hill=1.5;			// mean of Hill-coefficient (drawn randomly for each web from normal distribution)
double sigma_hill=0.2;		// sigma of distribution of Hill-coefficient
double hill;				// Hill coefficient of the while network (if applicable)

double x_resp=0.314;		// respiration rate
double C_interference;		// predator interference
double mu_Cint=1;			// mean predator interference
double sigma_Cint=0.1;		// sigma of distribution of predator interference
double lambda=0.85;			// assimilation efficiency

double B_b=10;				// intercept of initial basal species bodymass
double B_c=10;				// intercept of initial predator bodymass

//	parameters of nutrient model
double C_N[1];							// array for nutrient contents in plants (values are assigned in the main function)
double C1=1;							// content of first nutrient
double C2=0.5;							// content of second nutrient
double C3=0.25;							// content of all other nutrients
double D_N=0.25;						// nutrient turnover rate
double K_min=0.1;						// minimal nutrient uptake half saturation density
double K_rel=0.1;						// width of interval for nutrient uptake half saturation densities
double mu_S=10;//4;						// mean nutrient supply concentration
double sigma_S=2;						// standard deviation of nutrient supply concentration

//	default parameters for loops and flags
double def;					// default parameter
int GF=0;					// global flag, for testing purposes
int output_flag=0;			// 0 indicates statistical simulations, other number gives timeseries of a certain web
double cv_flag;				// indicates potentially pathological coefficients of variation


// functions:
double *web_calc(gsl_rng *r,FILE *data1);
int dynamics(double t, const double y[], double ydot[], void *params);
int pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mas);
int output(gsl_matrix *Ap,gsl_vector *mas,double B[],double iniB[],double meanB[],double VarCoeffB[],void *params,FILE *data1);
int Average_TL(gsl_matrix *Ap, double TLav[], double TLmin[]);
int Min_TL(gsl_matrix *Ap, double TLmin[]);
int Eff_TL(gsl_matrix *Ap, double meanB[], void *params, double TLeff[], double TLmin[]);
double functional_diversity(gsl_vector *mas);
int flows(double meanB[], void *params, double energyflows[]);
int mass_abundance(double meanB[], gsl_vector *mas, double fitresults[]);

//*******************************************************************
// Main Part
//*******************************************************************
int main(int argc, char* argv[]) // waits for input of one parameter (two-digit number to invoke simulation job)
{
	FILE *data1;

	int i,j;  // identifier of predator, i,  and prey, j 

	int zahl;
	char einer='0',zehner='0';
	const char * str="0";

	str=argv[1];							// reads input parameter (string, d.h. array von charactern)

	zehner = (char)str[0];					// gets first digit of input parameter (character)
	einer = (char)str[1];  				// gets second digit of input parameter (character)
	zahl = 10 * (zehner - '0') + (einer -'0');

	seed+=zahl;								// modifiziert die seed entsprechend des uebergebenen Parameters
	strcat(filename,str);					// haengt den uebergebenen Parameter an den Namen des Ausgabefiles an
	Web_Id=DivLev*Imax*(zahl-1);			// erhoeht die Web_Id

	gsl_rng_default_seed=seed;				// Startwert fuer den RNG
	gsl_rng *r=gsl_rng_alloc(gsl_rng_default);

	C_N[0]=C1;
	C_N[1]=C2;

	data1=fopen(filename,"w");				// creates empty file with the name given in the header of the code

	if(output_flag==0)						// header for webstats
	{
		fprintf(data1,"#job_id=%s,	seed=%d\n",str,seed);

		fprintf(data1,"WebID	S_c	hill	pred_int	P	P_b	P_c P_n	");
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
		fprintf(data1,"interc_b	expon_b	R2_b	interc_c	expon_c	R2_c	interc	expon	R2	");
		fprintf(data1,"cv_flag	\n"); 
	}
	else									// header for timeseries
	{
		fprintf(data1, "timestep ");
		for(i=0; i<S_b; i++) 
				fprintf(data1, "B%d ", i+1);
		for(i=0; i<S_c; i++) 
				fprintf(data1, "C%d ", i+1);
		for(i=0; i<N; i++) 
				fprintf(data1, "N%d ", i+1);	
		fprintf(data1, "\n");
	}
	fclose(data1);

	for(j=0;j<DivLev;j++)
	{
		if(output_flag==0)
			S_c=j*DivStep;

		S=S_c+S_b;

		for(i=0;i<Imax;i++)
		{
//			GF=i;
//			printf("%d	%d\n",j,i);
			Web_Id++;

			if(output_flag!=0)
			{
				if(Web_Id==output_flag)
					tend=30000;
				else
				{
					tend=0.1;
					teval=0.1;
				}
			}

			data1=fopen(filename,"a");		//re-opens output file to amend it
			web_calc(r,data1);
			fclose(data1);					// forces programm to write data to file
		}
	}

	gsl_rng_free(r);

	return(0);
}

//***************************************************************************************
// Funktion erzeugt ein ADBM-Netz und stellt die Umgebung fuer die Populationsdynamik
//*************************************************************************************** 
double *web_calc(gsl_rng *r,FILE *data1)
{
	int i,j,k;
	double temp1,temp2,temp3,temp4;
	gsl_matrix *Ap=gsl_matrix_calloc(S,S);						// adjacency matrix
	gsl_vector *mas=gsl_vector_calloc(S);						// mean bodymasses of species

	gsl_vector *defvec=gsl_vector_calloc(S);					// default vector: use it, but don't store important data in it!

	double *iniB=(double *) calloc(S+N,sizeof(double));			// initial biomass denisties
	double *meanB=(double *) calloc(S+N,sizeof(double));		// mean biomass densities 
	double *meanB2=(double *) calloc(S+N,sizeof(double));		// mean biomass densities (for comparison)
	double *VarB=(double *) calloc(S+N,sizeof(double));			// variances of biomasses
	double *VarCoeffB=(double *) calloc(S+N,sizeof(double));	// coefficients of variation

//	*********** generate network-topology  ***********
	pdef_structure(r,Ap,mas);									// calls function for network topology	

	if(output_flag==Web_Id)
	{
		for(i=0;i<S;i++)
		{
			for(j=0;j<S;j++)
				printf("%d ",(int)gsl_matrix_get(Ap,i,j));
			printf("	%lf\n",gsl_vector_get(mas,i));
		}
	}

//	*********** Dynamik berechnen: Parameter definieren, Differentialgleichung integrieren ************
	double *params=(double *) calloc(2*S*S+5*S+2*S_b*N+N+S_b,sizeof(double));		// parameters for population dynamics
	gsl_matrix *A=gsl_matrix_calloc(S,S);						// attack rates
	gsl_matrix *H=gsl_matrix_calloc(S,S);						// handling times

	gsl_matrix_memcpy(A,Ap);

 	double *TLmin=(double *) calloc(S,sizeof(double));
	Min_TL(Ap,TLmin);											// minimale trophische Level berechnen

	temp1=0;												
	while(temp1==0)
	{
		hill=gsl_ran_gaussian(r,sigma_hill)+mu_hill;			// Hill-Koeffizient konstant fuer das ganze Netzwerk
		if(1<=hill&&hill<=2.0)
			temp1=1;
	}

	temp1=0;												
	while(temp1==0)
	{
		C_interference=gsl_ran_gaussian(r,sigma_Cint)+mu_Cint;			// ziehung der Predatoren Interferenz konstant fuer das ganze Netzwerk
		if(0.75<=C_interference && C_interference<=1.25)
			temp1=1;
	}

	for(i=0;i<S;i++)
	{
		if(TLmin[i]!=1)
		{
			gsl_vector_view tempp=gsl_matrix_row(A,i);			// bug detected! hier wurde vorher Ap manipuliert!
			temp1=gsl_blas_dasum(&tempp.vector);
			gsl_vector_scale(&tempp.vector,1/temp1);

			temp1=gsl_vector_get(mas,i);

			for(j=0;j<S;j++)
			{
				temp2=gsl_vector_get(mas,j);
				temp3=pow(temp2,a_p)*pow(temp1,a_c)*pow((temp1/(temp2*R_opt))*exp(1-(temp1/(temp2*R_opt))),g);
				gsl_matrix_set(A,i,j,a_0*temp3*gsl_matrix_get(A,i,j));		// attack rates

				temp4=h_0*pow(temp1,h_c)*pow(temp2,h_p);
				gsl_matrix_set(H,i,j,temp4);								// handling times
			}
		}

		for(j=0;j<S;j++)
		{
			*(params+i*S+j)=gsl_matrix_get(A,i,j);
			*(params+S*S+i*S+j)=gsl_matrix_get(H,i,j);
		}

		*(params+2*S*S+i)=x_resp*pow(gsl_vector_get(mas,i),-0.25);			// Respirationsraten

		*(params+2*S*S+S+i)=C_interference;									// Interferenz-Konkurrenz

		if(TLmin[i]==1)
			*(params+2*S*S+2*S+i)=1;										// identifiziert Basalarten

/*		temp1=0;															// wenn diese Zeilen nicht auskommentiert sind, ist 
		while(temp1==0)														// der Hill-Koeffizient Beute-spezifisch
		{
			hill=gsl_ran_gaussian(r,sigma_hill)+mu_hill;					// Hill-Koeffizienten
			if(1<=hill&&hill<=2)
				temp1=1;
		}
*/		*(params+2*S*S+3*S+i)=hill;
		*(params+2*S*S+4*S+i)=gsl_vector_get(mas,i);						// Koerpermassen
	}
	for(i=0;i<S_b;i++)
	{
		for(j=0;j<N;j++)
		{
			*(params+2*S*S+5*S+i*N+j)=C_N[j];								// nutrient contents in plants
			*(params+2*S*S+5*S+S_b*N+i*N+j)=K_min+K_rel*gsl_rng_uniform(r);	// nutrient uptake half saturation densities
		}
	}
	for(i=0;i<N;i++)
	{
		temp1=gsl_ran_gaussian(r,sigma_S)+mu_S;
		if(temp1<=0)
			i--;
		else
			*(params+2*S*S+5*S+2*S_b*N+i)=temp1;							// nutrient supply concentrations
	}
	for(i=0;i<S_b;i++)
		*(params+2*S*S+5*S+2*S_b*N+N+i)=pow(gsl_vector_get(mas,i),-0.25);	// max. nutrient uptake rates

//	*********** initialise biomass densities and solve ODE *********** 
	double *B=(double *) calloc(S+N,sizeof(double));

	const gsl_odeiv_step_type *Solv=gsl_odeiv_step_rkf45;		// Algorithmus fuer DGL-Solver: Runge-Kutta-Fehlberg
	gsl_odeiv_control *c=gsl_odeiv_control_y_new(1e-8,1e-10);	// absolute und relative Fehlertoleranz
	gsl_odeiv_step *s=gsl_odeiv_step_alloc(Solv,S+N);
	gsl_odeiv_evolve *e=gsl_odeiv_evolve_alloc(S+N);	
	gsl_odeiv_system sys={dynamics,NULL,S+N,params};

	for(i=0;i<S;i++)
	{
		if(TLmin[i]==1)
			B[i]=B_b*gsl_rng_uniform_pos(r);					// Basalarten
		else
			B[i]=B_c*gsl_rng_uniform_pos(r);					// Konsumenten
	}
	for(i=0;i<N;i++)
	{
			temp1=(double) *(params+2*S*S+5*S+2*S_b*N+i);
			B[S+i]=0.5*temp1+0.5*temp1*gsl_rng_uniform_pos(r);	// nutrients
	}
	for(i=0;i<S+N;i++)
		iniB[i]=B[i];											// saving initial biomass values

	double t=0.0;												// Startzeit festlegen
	double h=1e-5;												// Anfangsschrittweite

	while(t<tend)
	{
		int status=gsl_odeiv_evolve_apply(e,c,s,&sys,&t,tend,&h,B);
		if(status!= GSL_SUCCESS)
			break;

		for(i=0;i<S+N;i++)
		{	
			if(B[i]<1e-6)										// Spezies mit zu kleiner Biomasse sterben aus 
				B[i]=0;	
		}
		
		if(output_flag==Web_Id)
		{
			fprintf(data1,"%lf ",t);
			for(i=0;i<S+N;i++)
				fprintf(data1,"%g ",B[i]);							// Zeitreihe der Biomassen ausgeben
			fprintf(data1,"\n");
		}
	}

	double tau=t,tau_old=0,Delta_tau;

	while(t<tend+teval)
	{
		int status=gsl_odeiv_evolve_apply(e,c,s,&sys,&t,tend+teval,&h,B);
		if(status!= GSL_SUCCESS)
			break;

		for(i=0;i<S;i++)
		{	
			if(B[i]<1e-6)										// Spezies mit zu kleiner Biomasse sterben aus 
				B[i]=0;	
		}

		tau_old=tau;
		tau=t;
		Delta_tau=tau-tau_old;

		for(i=0;i<S+N;i++)
			meanB[i]+=B[i]*Delta_tau;
	}
	for(i=0;i<S+N;i++)
		meanB[i]=meanB[i]/teval;								// Mittelwert der Biomassen aller Arten

	tau=t;
	gsl_vector *tempvec=gsl_vector_calloc(S+N); 
	gsl_vector_set_all(tempvec,1);
	
	while(t<tend+2*teval)
	{
		int status=gsl_odeiv_evolve_apply(e,c,s,&sys,&t,tend+2*teval,&h,B);
		if(status!= GSL_SUCCESS)
			break;

		for(i=0;i<S;i++)
		{	
			if(B[i]<1e-6)										// Spezies mit zu kleiner Biomasse sterben aus 
				B[i]=0;	
		}

		
		tau_old=tau;
		tau=t;
		Delta_tau=tau-tau_old;

		
		for(i=0;i<S+N;i++)
		{
			VarB[i]+=(B[i]-meanB[i])*(B[i]-meanB[i])*Delta_tau;
			meanB2[i]+=B[i]*Delta_tau;	
		}
					
	}
	temp1=0;
	temp2=0;
	for(i=0;i<S+N;i++)
	{
		meanB2[i]=meanB2[i]/teval;
		temp1+=meanB[i];
		temp2+=meanB2[i];
	}
	cv_flag=temp1-temp2;

	for(i=0;i<S+N;i++)
	{
		VarB[i]=VarB[i]/teval;		
		if(meanB[i]>0)
			VarCoeffB[i]=sqrt(VarB[i])/meanB[i];				// Coefficient of variation
		else
			VarCoeffB[i]=0;
	}

//	*********** Netzwerk auswerten und output generieren ************
	if(output_flag==0)
		output(Ap,mas,B,iniB,meanB,VarCoeffB,params,data1);

//	*********** Free memories ***********
	gsl_odeiv_step_free(s);
	gsl_odeiv_control_free(c);
	gsl_odeiv_evolve_free(e);
	gsl_matrix_free(Ap);
	gsl_matrix_free(A);
	gsl_matrix_free(H);
	gsl_vector_free(mas);
	gsl_vector_free(defvec);
	gsl_vector_free(tempvec);
	free(params);
	free(TLmin);
	free(meanB);
	free(meanB2);
	free(VarB);
	free(VarCoeffB);

	return 0;
}

//***************************************************************************************
//	Topologie eines koerpermassenbasierten Netzes erzeugen
//***************************************************************************************
int pdef_structure(gsl_rng *r,gsl_matrix *Ap, gsl_vector *mas)
{
	int i,j,flag=0;
	double Wkeit,sigma_i,zeta_act;
	double temp1,temp2,m_opt,a_max;

	while(flag==0)
	{
		flag=1;
		gsl_vector_set_all(mas,0);

//	********* Koerpermassen bestimmen *********
		zeta_act=zeta_c*log(10);			// Umrechnung auf ln des Mittelwertes
		for(i=0;i<S_c;i++)
		{
			temp1=gsl_ran_lognormal(r,zeta_act,sigma_c);	
			if(temp1>exp(zeta_act)/cutoff_c&&temp1<exp(zeta_act)*cutoff_c)		// minimale und maximale Masse
				gsl_vector_set(mas,S_b+i,temp1);								// adulte Masse
			else
				i--;
		}
		gsl_sort_vector(mas);

		gsl_vector_view mas_b_vec=gsl_vector_subvector(mas,0,S_b);
		gsl_vector *mas_b=&mas_b_vec.vector;

		zeta_act=zeta_b*log(10);			// Umrechnung auf ln des Mittelwertes
		for(i=0;i<S_b;i++)
		{
			temp1=gsl_ran_lognormal(r,zeta_act,sigma_b);	
			if(temp1>exp(zeta_act)/cutoff_b&&temp1<exp(zeta_act)*cutoff_b)		// minimale und maximale Masse
				gsl_vector_set(mas_b,i,temp1);								// adulte Masse
			else
				i--;
		}
		gsl_sort_vector(mas_b);

//	********** Matrix Ap befuellen: Rickers attack rate *************
		for(i=0;i<S_c;i++)
		{
			temp1=gsl_vector_get(mas,S_b+i);
			a_max=pow(gsl_vector_get(mas,S_b+i)/R_opt,0.25);

			for(j=0;j<S;j++)
			{
				temp2=gsl_vector_get(mas,j);
				if(pow(temp2,0.25)*pow((temp1/(temp2*R_opt))*exp(1-(temp1/(temp2*R_opt))),g)>=cutoff2*a_max)
					gsl_matrix_set(Ap,S_b+i,j,1);
			}

			gsl_vector_view tempp=gsl_matrix_row(Ap,S_b+i);
			flag=flag*(1-gsl_vector_isnull(&tempp.vector));
		}

	}

	return(0);
}


//***************************************************************************************
//	Output generieren
//***************************************************************************************
int output(gsl_matrix *Ap, gsl_vector *mas, double B[], double iniB[], double meanB[], double VarCoeffB[], void *params, FILE *data1)
{
	int i,j,k;
	gsl_vector *defvec=gsl_vector_calloc(S);				// default vector: use it, but don't store important data in it!
	gsl_matrix *defmat=gsl_matrix_calloc(S,S);				// default matrix: same as with defvec....
	double *pparams=(double *)params;						// pointer to first Element of params...

//	Persistenz, mittlere Biomassen, Variationskoeffizienten

	double persistence=0,persistence_b=0,persistence_c=0,persistence_n=0;	// Persistenz (total, basal, Konsumer)
	double meanB_n=0,meanB_b=0,meanB_c=0,meanB_total=0;			// mittlere Biomassendichten
	double CV_n=0,CV_b=0,CV_c=0,CV=0;										// Variationskoeffizienten

	for(i=0;i<S;i++)
	{
		if(B[i]>=1e-6)
		{
			persistence++;
			meanB_total+=meanB[i];
			CV+=VarCoeffB[i];
		}
	}
	for(i=0;i<S_b;i++)
	{
		if(B[i]>=1e-6)
		{
			persistence_b++;
			meanB_b+=meanB[i];
			CV_b+=VarCoeffB[i];
		}
	}
	for(i=0;i<S_c;i++)
	{
		if(B[S_b+i]>=1e-6)
		{
			persistence_c++;
			meanB_c+=meanB[S_b+i];
			CV_c+=VarCoeffB[S_b+i];
		}
	}
	for(i=0;i<N;i++)
	{
		if(B[S+i]>=1e-6)
		{
			persistence_n++;
			meanB_n+=meanB[S+i];
			CV_n+=VarCoeffB[S+i];
		}
	}

	if(persistence!=0)
	{
		meanB_total=meanB_total/persistence;
		CV=CV/persistence;
	}
	else
	{
		meanB_total=0;
		CV=0;
	}
	if(persistence_b!=0)
	{
		meanB_b=meanB_b/persistence_b;
		CV_b=CV_b/persistence_b;
	}
	else
	{
		meanB_b=0;
		CV_b=0;
	}
	if(persistence_c!=0)
	{
		meanB_c=meanB_c/persistence_c;
		CV_c=CV_c/persistence_c;
	}
	else
	{
		meanB_c=0;
		CV_c=0;
	}
	if(persistence_n!=0)
	{
		meanB_n=meanB_n/persistence_n;
		CV_n=CV_n/persistence_n;
	}
	else
	{
		meanB_n=0;
		CV_n=0;
	}

	persistence=persistence/S;
	persistence_b=persistence_b/S_b;
	persistence_c=persistence_c/S_c;
	persistence_n=persistence_n/N;

//	fuer finale Groessen: ausgestorbene Spezies aus dem Netzwerk entfernen

	gsl_matrix_memcpy(defmat,Ap);
	gsl_vector_memcpy(defvec,mas);
	
	for(i=0;i<S;i++)
	{
		if(B[i]<1e-6)
		{
			gsl_vector_set(defvec,i,0);
			for(j=0;j<S;j++)
			{
				gsl_matrix_set(defmat,i,j,0);
				gsl_matrix_set(defmat,j,i,0);
			}
		}
	}

//	trophische Level berechnen

 	double *TLmin=(double *) calloc(S,sizeof(double));		// Speicher fuer troph. Level (Def: distance)
	double *TLav=(double *) calloc(S,sizeof(double));		// Speicher fuer troph. Level (Def: gemittelt)
	double *TLeff=(double *) calloc(S,sizeof(double));		// Speicher fuer troph. Level (Def: effective)

	double meanTLmin_ini=0,meanTLmin_fin=0;
	double meanTLav_ini=0,meanTLav_fin=0;
	double meanTLeff=0;
	int iniOmn=0, finOmn=0;
	
	Min_TL(Ap,TLmin);										// initiale minimale trophische Level berechnen
	Average_TL(Ap,TLav,TLmin);								// initiale mittlere trophische Level berechnen

	for(i=0;i<S_c;i++)
	{
		meanTLmin_ini+=TLmin[S_b+i];
		meanTLav_ini+=TLav[S_b+i];
	}
	meanTLmin_ini=meanTLmin_ini/S_c;						// mittleres initiales minimales trophisches Level der Konsumenten
	meanTLav_ini=meanTLav_ini/S_c;							// mittleres initiales gemitteltes trophisches Level der Konsumenten


	for(i = 0; i < S_c; i++) 
	{
		k=0;
		for(j = 0; j < S; j++) 
		{		
			if(gsl_matrix_get(Ap,S_b+i,j) == 1 && TLmin[j] >= TLmin[S_b+i])
				k=1;										
		}
		iniOmn+=k;										// initial number of omnivores
	}

	Min_TL(defmat,TLmin);									// finale minimale trophische Level berechnen
	Average_TL(defmat,TLav,TLmin);							// finale mittlere trophische Level berechnen
	Eff_TL(defmat,meanB,params,TLeff,TLmin);				// effective trophische Level berechnen

	for(i=0;i<S_c;i++)
	{
		if(B[S_b+i]>=1e-6)
		{
			meanTLmin_fin+=TLmin[S_b+i];
			meanTLav_fin+=TLav[S_b+i];
			meanTLeff+=TLeff[S_b+i];
		}
	}
	if(persistence_c!=0)
	{
		meanTLmin_fin=meanTLmin_fin/(S_c*persistence_c);	// mittleres finales minimales trophisches Level der Konsumenten
		meanTLav_fin=meanTLav_fin/(S_c*persistence_c);		// mittleres finales gemitteltes trophisches Level der Konsumenten
		meanTLeff=meanTLeff/(S_c*persistence_c);			// mittleres effektives trophischesl Level der Konsumenten
	}
	else
	{
		meanTLmin_fin=0;
		meanTLav_fin=0;
		meanTLeff=0;
	}
	
	for(i = 0; i < S_c; i++) 
	{
		k = 0;
		for(j = 0; j < S; j++) 
		{		
			if(gsl_matrix_get(defmat,S_b+i,j) == 1 && TLmin[j] >= TLmin[S_b+i]) 
				k = 1;
		}
		finOmn += k;										// final number of omnivores	
	}
		
//	Links,  Massen-, Vulnerabilitaets- und Generalitaetsverteilungen

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
	
	for(i=0; i<S; i++)
	{
		for(j=0; j<S; j++)
		{
			ini_l+=gsl_matrix_get(Ap,i,j);									//initial Number of Links
			fin_l+=gsl_matrix_get(defmat,i,j);								//final Number of Links
		}
	}

	for(i=0;i<S;i++)
	{
		logmas[i]=log10(gsl_vector_get(mas,i));
	
		gsl_vector_view tempp_in=gsl_matrix_row(Ap,i);
		gen[i]=gsl_blas_dasum(&tempp_in.vector);							// generality Verteilung

		gsl_vector_view tempp_out=gsl_matrix_column(Ap,i);
		vul[i]=gsl_blas_dasum(&tempp_out.vector);							// out-degree Verteilung
	}
	
	for(i=0;i<S_b;i++)
	{
		inivul_b+=vul[i]/S_b;												// mittlere initiale Vulnerabilitaet der Basalarten
		inilogmas_b+=logmas[i]/S_b;											// mittlere initiale log. Koerpermasse der Basalarten
		inidens_b+= iniB[i]/gsl_vector_get(mas,i);							// initiale Gesamt-Individuendichte aller Basalarten 
	}
	for(i=0;i<S_c;i++)
	{
		inivul_c+=vul[S_b+i]/S_c;											// mittlere initiale Vulnerabilitaet der Konsumenten
		inigen_c+=gen[S_b+i]/S_c;											// mittlere initiale Generalitaet der Konsumenten
		inilogmas_c+=logmas[S_b+i]/S_c;										// mittlere initiale log. Koerpermasse der Konsumenten
		inidens_c+= iniB[S_b+i]/gsl_vector_get(mas,S_b+i);					// initiale Gesamt-Individuendichte aller Konsumenten
	}
	
	 for(i=0;i<S_c;i++)
	{
			if(vul[S_b+i] == 0 && gen[S_b+i] > 0) 
			initop_c+=1;													// initial Number of Top predators
	}
	
	for(i=0;i<S_b;i++)
	{
			if(vul[i] == 0) 
			initop_b+=1;													// initial Number of predator-free basal species 
	}
	
	for(i=0;i<S_b;i++)
	{
		iniSDvul_b+=(vul[i]-inivul_b)*(vul[i]-inivul_b)/S_b;
		iniSDlogmas_b+=(logmas[i]-inilogmas_b)*(logmas[i]-inilogmas_b)/S_b;
	}
	for(i=0;i<S_c;i++)
	{
		iniSDvul_c+=(vul[S_b+i]-inivul_c)*(vul[S_b+i]-inivul_c)/S_c;
		iniSDgen_c+=(gen[S_b+i]-inigen_c)*(gen[S_b+i]-inigen_c)/S_c;
		iniSDlogmas_c+=(logmas[S_b+i]-inilogmas_c)*(logmas[S_b+i]-inilogmas_c)/S_c;
	}
	iniSDvul_b=sqrt(iniSDvul_b);											// SD initiale Vulnerabilitaet der Basalarten
	iniSDlogmas_b=sqrt(iniSDlogmas_b);										// SD initiale log. Koerpermasse der Basalarten
	iniSDvul_c=sqrt(iniSDvul_c);											// SD initiale Vulnerabilitaet der Konsumenten
	iniSDgen_c=sqrt(iniSDgen_c);											// SD initiale Generalitaet der Konsumenten
	iniSDlogmas_c=sqrt(iniSDlogmas_c);										// SD initiale log. Koerpemasse der Konsumenten

//	finale Groessen:

	for(i=0;i<S;i++)
	{
		gsl_vector_view tempp_in=gsl_matrix_row(defmat,i);
		gen[i]=gsl_blas_dasum(&tempp_in.vector);							// generality Verteilung

		gsl_vector_view tempp_out=gsl_matrix_column(defmat,i);
		vul[i]=gsl_blas_dasum(&tempp_out.vector);							// out-degree Verteilung
	}

	for(i=0;i<S_b;i++)
	{
		if(B[i]>=1e-6)
		{
			finvul_b+=vul[i]/(S_b*persistence_b);							// mittlere finale Vulnerabilitaet der Basalarten
			finlogmas_b+=logmas[i]/(S_b*persistence_b);						// mittlere finale log. Koerpermasse der Basalarten
			findens_b+= meanB[i]/gsl_vector_get(mas,i);						// finale Gesamt-Individuendichte aller Basalarten
		}
	}
	for(i=0;i<S_c;i++)
	{
		if(B[S_b+i]>=1e-6)
		{
			finvul_c+=vul[S_b+i]/(S_c*persistence_c);						// mittlere finale Vulnerabilitaet der Konsumenten
			fingen_c+=gen[S_b+i]/(S_c*persistence_c);						// mittlere finale Generalitaet der Konsumenten
			finlogmas_c+=logmas[S_b+i]/(S_c*persistence_c);					// mittlere finale log. Koerpermasse der Konsumenten
			findens_c+= meanB[S_b+i]/gsl_vector_get(mas,S_b+i);				// finale Gesamt-Individuendichte aller Konsumenten
		}
	}
	
	for(i=0;i<S_c;i++)
	{
		if(B[i]>=1e-6)
			if(vul[S_b+i] == 0 && gen[S_b+i] > 0) 
				fintop_c+=1;												// final Number of Top predators
	}
	
	for(i=0;i<S_b;i++)
	{
		if(B[i]>=1e-6)
			if(vul[i] == 0) 
				fintop_b+=1;												// final Number of predator-free basal species 
	}
	
	for(i=0;i<S_b;i++)
	{
		if(B[i]>=1e-6)
		{
			finSDvul_b+=(vul[i]-finvul_b)*(vul[i]-finvul_b)/(S_b*persistence_b);
			finSDlogmas_b+=(logmas[i]-finlogmas_b)*(logmas[i]-finlogmas_b)/(S_b*persistence_b);
		}
	}
	for(i=0;i<S_c;i++)
	{
		if(B[S_b+i]>=1e-6)
		{
			finSDvul_c+=(vul[S_b+i]-finvul_c)*(vul[S_b+i]-finvul_c)/(S_c*persistence_c);
			finSDgen_c+=(gen[S_b+i]-fingen_c)*(gen[S_b+i]-fingen_c)/(S_c*persistence_c);
			finSDlogmas_c+=(logmas[S_b+i]-finlogmas_c)*(logmas[S_b+i]-finlogmas_c)/(S_c*persistence_c);
		}
	}
	finSDvul_b=sqrt(finSDvul_b);											// SD finale Vulnerabilitaet der Basalarten
	finSDlogmas_b=sqrt(finSDlogmas_b);										// SD finale log. Koerpermasse der Basalarten
	finSDvul_c=sqrt(finSDvul_c);											// SD finale Vulnerabilitaet der Konsumenten
	finSDgen_c=sqrt(finSDgen_c);											// SD finale Generalitaet der Konsumenten
	finSDlogmas_c=sqrt(finSDlogmas_c);										// SD finale log. Koerpemasse der Konsumenten


//	funktionelle Diversitaet

	double initial_funct_div,final_funct_div;					// initiale und finale funktionelle Diversitaet

	initial_funct_div=functional_diversity(mas);
	final_funct_div=functional_diversity(defvec);

//	Energiefluesse

	double *energyflows=(double *) calloc(6,sizeof(double));	// 0: nutrient influx, 1: nutrient-basal, 2: basal-consumer, 3: igp, 4: basal resp., 5: cons. resp.

	flows(meanB,params,energyflows);

//	mass-abundance relation

	double *fitresults=(double *) calloc(9,sizeof(double));

	mass_abundance(meanB,mas,fitresults);

//	Ausgabe in Datei

	fprintf(data1,"%d	%d	%.6g	%.6g	%.6g	%.6g	%.6g	%.6g	",Web_Id, S_c, hill, C_interference, persistence, persistence_b, persistence_c, persistence_n);
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
	for(i=0;i<9;i++)
		fprintf(data1,"%.6g	",fitresults[i]);
	fprintf(data1,"%.6g\n", cv_flag);


	free(TLmin);
	free(TLav);
	free(TLeff);
	free(vul);
	free(gen);
	free(fitresults);
	gsl_vector_free(defvec);
	gsl_matrix_free(defmat);

	return(0);
}


//***************************************************************************************
// Energiefluesse berechnen
//***************************************************************************************
int flows(double meanB[], void *params, double energyflows[])
{
	int i,j;
	double *pparams=(double *)params;									// Pointer zum ersten Element von params

	gsl_vector *TBvec=gsl_vector_calloc(S+N);							// all biomasses and nutrient concentrations	
	gsl_vector_const_view TB_vec1=gsl_vector_const_view_array(meanB,S+N);
	gsl_vector_memcpy(TBvec,&TB_vec1.vector);

	gsl_vector_view B_vec=gsl_vector_subvector(TBvec,0,S);				// all biomasses (plants + animals)
	gsl_vector *Bvec=&B_vec.vector;

	gsl_vector_view PB_vec=gsl_vector_subvector(TBvec,0,S_b);			// Basal biomasses
	gsl_vector *PBvec=&PB_vec.vector;

	gsl_vector_view NB_vec=gsl_vector_subvector(TBvec,S,N);				// nutrient concentrations
	gsl_vector *NBvec=&NB_vec.vector;

	gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);				// a_ij*f_ij
	gsl_matrix *Amat=&A_mat.matrix;

	gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);			// H_ij
	gsl_matrix *Hmat=&H_mat.matrix;

	gsl_vector_view R_vec=gsl_vector_view_array(pparams+2*S*S,S);			// respiration rates
	gsl_vector *Rvec=&R_vec.vector;

	gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);		// interference vector
	gsl_vector *Cvec=&C_vec.vector;

	gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);		// hill coefficients
	gsl_vector *Hvec=&H_vec.vector;	

	gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);		// body masses
	gsl_vector *Mvec=&M_vec.vector;	

	gsl_matrix_view C_mat=gsl_matrix_view_array(pparams+2*S*S+5*S,S_b,N);	// nutrient contents in plants
	gsl_matrix *Cmat=&C_mat.matrix;											// (P rows, N columns)

	gsl_matrix_view K_mat=gsl_matrix_view_array(pparams+2*S*S+5*S+S_b*N,S_b,N);	// nutrient uptake half saturation densities
	gsl_matrix *Kmat=&K_mat.matrix;											// (P rows, N columns)

	gsl_vector_view S_vec=gsl_vector_view_array(pparams+2*S*S+5*S+2*S_b*N,N);	// nutrient supply concentrations
	gsl_vector *Svec=&S_vec.vector;	

	gsl_vector_view U_vec=gsl_vector_view_array(pparams+2*S*S+5*S+2*S_b*N+N,S_b);// max. nutrient uptake rates
	gsl_vector *Uvec=&U_vec.vector;	

	gsl_vector *Prey=gsl_vector_calloc(S);									// prey biomasses
	gsl_vector *Dvec=gsl_vector_calloc(S);									// death (mortality+predation losses)

	gsl_matrix *mmat=gsl_matrix_calloc(S,S);
	gsl_vector *rvec=gsl_vector_calloc(S);
	gsl_vector *svec=gsl_vector_calloc(S);
	gsl_vector *nvec=gsl_vector_calloc(N);
	gsl_vector *mvec=gsl_vector_calloc(N);

//	************ in case of Type 3 functional response ***********
	gsl_vector_memcpy(Prey,Bvec);

	for(i=0;i<S;i++)
	{
		if(gsl_vector_get(Prey,i)<0)
			gsl_vector_set(Prey,i,0);				// paranoia (this doesn't mean that you can safely remove these lines!)

		gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));
	}

//	************* functional responses zusammenbauen *************
	gsl_matrix_memcpy(mmat,Hmat);
	gsl_matrix_mul_elements(mmat,Amat);

	gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);		// handling time term: r_i=sum_j a_ij*h_ij*Prey_j
	gsl_vector_memcpy(svec,Bvec);
	gsl_vector_mul(svec,Cvec);								// predator interference: s_i=C_i*B_i
	gsl_vector_add(rvec,svec);
	gsl_vector_add_constant(rvec,1);
	gsl_vector_mul(rvec,Mvec);								// rvec: denominator of functional response

	gsl_vector_memcpy(svec,Bvec);							// predator biomass
	gsl_vector_div(svec,rvec);								// ... mit Nenner der functional response multipliziert
	gsl_blas_dgemv(CblasTrans,1,Amat,svec,0,Dvec);			// Praedationsverluste pro Individuum:
	gsl_vector_mul(Dvec,Prey);								// multiply with target species' biomass

//	************* Terme fuer Basalarten zufuegen *************
	gsl_vector_view u_vec=gsl_vector_subvector(rvec,0,S_b);		// basal species growth rates
	gsl_vector *uvec=&u_vec.vector;
	for(i=0;i<S_b;i++)
	{
		gsl_vector_view K_vec=gsl_matrix_row(Kmat,i);
		gsl_vector *Kvec=&K_vec.vector;							// Kvec is a vector of length N (Kvec_j=Kmat_ij)

		gsl_vector_memcpy(nvec,NBvec);
		gsl_vector_add(nvec,Kvec);
		gsl_vector_memcpy(mvec,NBvec);
		gsl_vector_div(mvec,nvec);								// m_j=NB_j/(Kvec_j+NB_j)

		gsl_vector_set(uvec,i,gsl_vector_min(mvec));			// remember: uvec is just a vector_view of rvec!
	}
	gsl_vector_mul(uvec,Uvec);
	gsl_vector_mul(uvec,PBvec);									// uvec_i=PB_i*mas_i^-0.25*Min_j(NB_j/(NB_j+Kmat_ij))

//	************ flux 0: nutrient influx ***********
	energyflows[0]=D_N*gsl_blas_dasum(Svec);

//	************ flux 1: nutrients to basal species ***********
	energyflows[1]=gsl_blas_dasum(uvec);

//	************ flux 2: basal to consumers ***********
	for(i=0;i<S_b;i++)
		energyflows[2]+=gsl_vector_get(Dvec,i);

//	************ flux 3: consumers intraguild ***********
	for(i=0;i<S_c;i++)
		energyflows[3]+=gsl_vector_get(Dvec,S_b+i);

//	************ flux 4: basal respiration ***********
	for(i=0;i<S_b;i++)
		energyflows[4]+=meanB[i]*gsl_vector_get(Rvec,i);

//	************ flux 5: consumer respiration ***********
	for(i=0;i<S_c;i++)
		energyflows[5]+=meanB[S_b+i]*gsl_vector_get(Rvec,S_b+i);

//	********** free memory ************
	gsl_vector_free(TBvec);
	gsl_vector_free(Prey);
	gsl_vector_free(Dvec);
	gsl_vector_free(svec);
	gsl_vector_free(rvec);
	gsl_vector_free(nvec);
	gsl_vector_free(mvec);
	gsl_matrix_free(mmat);

	return(0);
}


//***************************************************************************************
// minimalen Abstand der Spezies von den Ressourcen bestimmen 
//***************************************************************************************
int Min_TL(gsl_matrix *Ap,double TLmin[])
{
	int i,j,k,counter;

	for(i=0;i<S;i++)
	{
		gsl_vector_view tempp=gsl_matrix_row(Ap,i);	
		if(gsl_vector_isnull(&tempp.vector)==1)
			TLmin[i]=1;											// Spezies auf dem ersten Level bestimmen
		else
			TLmin[i]=0;
	}	
	k=2;
	while(k<=S)			
	{
		counter=0;
		for(i=0;i<=S-1;i++)
		{
			if(TLmin[i]==0)
			{
				for(j=0;j<=S-1;j++)
				{
					if(TLmin[j]==k-1&&gsl_matrix_get(Ap,i,j)!=0)
					{	
						TLmin[i]=k;								// Spezies auf Level k, k>=2, bestimmen
						counter=counter+1;
					}
				}
			}
		}
		if(counter==0)
			k=S;	
			k=k+1;
	}  

	return(0);
}

//***************************************************************************************
//	mittleren Abstand der Spezies von den Ressourcen bestimmen
//***************************************************************************************
int Average_TL(gsl_matrix *Ap, double TLav[], double TLmin[])
{
// Loese Gleichung der Struktur (Matrix)Ta*(Vektor)xa=(Vektor)b 
	int i,j;
	double norm;
	int sig=1;
	gsl_vector *b=gsl_vector_calloc(S);					// Inhomogenitaet 
	gsl_vector *xa=gsl_vector_calloc(S);				// gesuchte Variablen
	gsl_permutation *p=gsl_permutation_calloc(S);	// Permutationsvektor
	gsl_matrix *Ta=gsl_matrix_calloc(S,S);				// Matrix des Gleichungssystems

	gsl_matrix_memcpy(Ta,Ap);

	for(i=0;i<S;i++)
	{
		gsl_vector_view tempp=gsl_matrix_row(Ta,i);
		norm=gsl_blas_dasum(&tempp.vector);
		
		if(norm!=0)
		{
			for(j=0;j<S;j++)
				gsl_matrix_set(Ta,i,j,gsl_matrix_get(Ta,i,j)/norm);
		}

		if(TLmin[i]==0)
		{
			for(j=0;j<S;j++)	
				gsl_matrix_set(Ta,i,j,0);
		}

	}
	for(i=0;i<S;i++)
		gsl_matrix_set(Ta,i,i,gsl_matrix_get(Ta,i,i)-1);

	for(i=0;i<S;i++)
	{
		if(TLmin[i]!=0)
			gsl_vector_set(b,i,-1);
	}

	gsl_linalg_LU_decomp(Ta,p,&sig);						// Matrix Ta zerlegen... 
	gsl_linalg_LU_solve(Ta,p,b,xa);						// ...und System loesen

	for(i=0;i<S;i++)
		TLav[i]=gsl_vector_get(xa,i);

	gsl_matrix_free(Ta);
	gsl_vector_free(xa);
	gsl_permutation_free(p);
	gsl_vector_free(b);

	return(0);
}


//***************************************************************************************
//	effektiven Abstand der Spezies von den Ressourcen bestimmen
//***************************************************************************************
int Eff_TL(gsl_matrix *Ap, double meanB[], void *params, double TLeff[], double TLmin[])
{
	int i,j;
	double *pparams=(double *)params;									// Pointer zum ersten Element von params

	gsl_vector *Bvec=gsl_vector_calloc(S);								// mean biomasses	
	gsl_vector_const_view B_vec1=gsl_vector_const_view_array(meanB,S);
	gsl_vector_memcpy(Bvec,&B_vec1.vector);

	gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);			// a_ij*f_ij
	gsl_matrix *Amat=&A_mat.matrix;

	gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);		// H_ij
	gsl_matrix *Hmat=&H_mat.matrix;

	gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);	// interference vector
	gsl_vector *Cvec=&C_vec.vector;

	gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);	// hill coefficients
	gsl_vector *Hvec=&H_vec.vector;	

	gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);	// body masses
	gsl_vector *Mvec=&M_vec.vector;	

	gsl_vector *Prey=gsl_vector_calloc(S);								// prey biomasses

	gsl_matrix *Effmat=gsl_matrix_calloc(S,S);
	gsl_matrix *mmat=gsl_matrix_calloc(S,S);
	gsl_vector *rvec=gsl_vector_calloc(S);
	gsl_vector *svec=gsl_vector_calloc(S);

//	************ in case of Type 3 functional response ***********
	gsl_vector_memcpy(Prey,Bvec);

	for(i=0;i<S;i++)
	{
		if(gsl_vector_get(Prey,i)<0)
			gsl_vector_set(Prey,i,0);				// paranoia (this doesn't mean that you can safely remove these lines!)

		gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));
	}

//	************* functional responses zusammenbauen *************
	gsl_matrix_memcpy(mmat,Hmat);
	gsl_matrix_mul_elements(mmat,Amat);

	gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);		// handling time term: r_i=sum_j a_ij*h_ij*Prey_j
	gsl_vector_memcpy(svec,Bvec);
	gsl_vector_mul(svec,Cvec);								// predator interference: s_i=C_i*B_i
	gsl_vector_add(rvec,svec);
	gsl_vector_add_constant(rvec,1);
	gsl_vector_mul(rvec,Mvec);								// rvec: denominator of functional response

	for(i=0;i<S;i++)
	{
		for(j=0;j<S;j++)
			gsl_matrix_set(Effmat,i,j,gsl_matrix_get(Amat,i,j)*gsl_vector_get(Prey,j)/gsl_vector_get(rvec,i));
	}

	Average_TL(Effmat, TLeff, TLmin);

//	********** free memory ************
	gsl_vector_free(Bvec);
	gsl_vector_free(Prey);
	gsl_vector_free(svec);
	gsl_vector_free(rvec);
	gsl_matrix_free(mmat);
	gsl_matrix_free(Effmat);

	return(0);
}

//***************************************************************************************
//	determines functional diversity
//***************************************************************************************
double functional_diversity(gsl_vector *mas)
{
	int i,j;
	double max;				// Maximum der Attackraten
	double M;				// Koerpermasse einer Art
	double M_prey;			// Koerpermasse der Beutespezies
	double attackrate;		// Attackrate einer Art 

	double L=12;			// log_10 Laenge der Koerpermassenachse
	double N=1200;			// Anzahl der Stuetzstellen, die ausgewertet werden
	double l=L/N;			// Laenge der Einzelintervalle
	double s;				// Stuetzstelle, an der Funktionswerte berechnet werden
	double I=0;				// Wert des Integrals (Mass fuer die funktionelle Diversitaet)

	for(j=0;j<N;j++)
	{		
		max=0;
		s=(j+0.5)*l;
		M_prey=pow(10,s);

		for(i=0;i<S_c;i++)
		{
			M=gsl_vector_get(mas,i+S_b);
			attackrate=pow(M/(M_prey*R_opt)*exp(1-(M/(M_prey*R_opt))),g);

			if(attackrate<cutoff2)
				attackrate=0;

			if(attackrate>max)
				max=attackrate;
		}

		I+=max*l;
	}
	I=I/L;					// Normierung

	return(I);
}


//***************************************************************************************
//	determine mass-abundance relation
//***************************************************************************************
int mass_abundance(double meanB[], gsl_vector *mas, double fitresults[])
{
	int i,j;
	int c=0,b=0;

	double exponent,exponent_b,exponent_c,intercept,intercept_b,intercept_c,R2,R2_b,R2_c;

	for(i=0;i<S_b;i++)
	{
		if(meanB[i]>0)
			b++;
	}
	for(i=0;i<S_c;i++)
	{
		if(meanB[i+S_b]>0)
			c++;
	}

	double *logmas_b=(double *) calloc(b,sizeof(double));
	double *logmas_c=(double *) calloc(c,sizeof(double));
	double *logmas=(double *) calloc(b+c,sizeof(double));
	double *logabundance_b=(double *) calloc(b,sizeof(double));
	double *logabundance_c=(double *) calloc(c,sizeof(double));
	double *logabundance=(double *) calloc(b+c,sizeof(double));

//	Daten: Logarithmen der Koerpermassen und Biomassedichten der Basalarten und der Konsumenten
	j=0;
	for(i=0;i<S_b;i++)
	{
		if(meanB[i]>0)
		{
			logmas_b[j]=log(gsl_vector_get(mas,i));			// Basalarten
			logabundance_b[j]=log(meanB[i]);
			j++;
		}
	}
	j=0;
	for(i=0;i<S_c;i++)
	{
		if(meanB[i+S_b]>0)
		{
			logmas_c[j]=log(gsl_vector_get(mas,i+S_b));		// Konsumenten
			logabundance_c[j]=log(meanB[i+S_b]);
			j++;
		}
	}
	for(i=0;i<b;i++)
	{
		logmas[i]=logmas_b[i];
		logabundance[i]=logabundance_b[i];					// Basalarten und Konsumenten zusammen
	}
	for(i=0;i<c;i++)
	{
		logmas[b+i]=logmas_c[i];
		logabundance[b+i]=logabundance_c[i];
	}

//	Fit eines linearen Modells an die Daten
	double meanX=0,meanY=0,SSXY=0,SSX=0,SS_tot=0,SS_err=0;
	for(i=0;i<b;i++)										// zunaechst die Basalarten...
	{
		meanX+=logmas_b[i]/b;
		meanY+=logabundance_b[i]/b;
	}
	for(i=0;i<b;i++)
	{
		SSXY+=(logmas_b[i]-meanX)*(logabundance_b[i]-meanY);
		SSX+=(logmas_b[i]-meanX)*(logmas_b[i]-meanX);
	}
	exponent_b=SSXY/SSX;
	intercept_b=meanY-exponent_b*meanX;

	for(i=0;i<b;i++)										// Berechnung des R^2-Wertes
	{
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

	for(i=0;i<c;i++)										// dann die Konsumenten...
	{
		meanX+=logmas_c[i]/c;
		meanY+=logabundance_c[i]/c;
	}
	for(i=0;i<c;i++)
	{
		SSXY+=(logmas_c[i]-meanX)*(logabundance_c[i]-meanY);
		SSX+=(logmas_c[i]-meanX)*(logmas_c[i]-meanX);
	}
	exponent_c=SSXY/SSX;
	intercept_c=meanY-exponent_c*meanX;

	for(i=0;i<c;i++)
	{
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

	for(i=0;i<b+c;i++)										// und dann alle Arten zusammen
	{
		meanX+=logmas[i]/(b+c);
		meanY+=logabundance[i]/(b+c);
	}
	for(i=0;i<b+c;i++)
	{
		SSXY+=(logmas[i]-meanX)*(logabundance[i]-meanY);
		SSX+=(logmas[i]-meanX)*(logmas[i]-meanX);
	}
	exponent=SSXY/SSX;
	intercept=meanY-exponent*meanX;

	for(i=0;i<b+c;i++)
	{
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
	
	return(0);
}


//***************************************************************************************
//	defines population dynamics equations
//***************************************************************************************
int dynamics(double t, const double B[], double Bdot[], void *params)
{
	int i,j;
	double *pparams=(double *)params;									// Pointer zum ersten Element von params

	gsl_vector *TBvec=gsl_vector_calloc(S+N);							// all biomasses and nutrient concentrations	
	gsl_vector_const_view TB_vec1=gsl_vector_const_view_array(B,S+N);
	gsl_vector_memcpy(TBvec,&TB_vec1.vector);

	gsl_vector_view B_vec=gsl_vector_subvector(TBvec,0,S);				// all biomasses (plants + animals)
	gsl_vector *Bvec=&B_vec.vector;

	gsl_vector_view PB_vec=gsl_vector_subvector(TBvec,0,S_b);			// Basal biomasses
	gsl_vector *PBvec=&PB_vec.vector;

	gsl_vector_view NB_vec=gsl_vector_subvector(TBvec,S,N);				// nutrient concentrations
	gsl_vector *NBvec=&NB_vec.vector;

	gsl_vector_view Bdot_vec=gsl_vector_view_array(Bdot,S);				// Vektor fuer Zeitableitungen
	gsl_vector *Bdotvec=&Bdot_vec.vector;

	gsl_vector_view NBdot_vec=gsl_vector_view_array(Bdot+S,N);				// Vektor fuer Zeitableitungen
	gsl_vector *NBdotvec=&NBdot_vec.vector;

	gsl_matrix_view A_mat=gsl_matrix_view_array(pparams,S,S);				// a_ij*f_ij
	gsl_matrix *Amat=&A_mat.matrix;

	gsl_matrix_view H_mat=gsl_matrix_view_array(pparams+S*S,S,S);			// H_ij
	gsl_matrix *Hmat=&H_mat.matrix;

	gsl_vector_view R_vec=gsl_vector_view_array(pparams+2*S*S,S);			// respiration rates
	gsl_vector *Rvec=&R_vec.vector;

	gsl_vector_view C_vec=gsl_vector_view_array(pparams+2*S*S+1*S,S);		// interference vector
	gsl_vector *Cvec=&C_vec.vector;

	gsl_vector_view P_vec=gsl_vector_view_array(pparams+2*S*S+2*S,S);		// basal vector (1 for basal species, 0 otherwise)
	gsl_vector *Pvec=&P_vec.vector;

	gsl_vector_view H_vec=gsl_vector_view_array(pparams+2*S*S+3*S,S);		// hill coefficients
	gsl_vector *Hvec=&H_vec.vector;	

	gsl_vector_view M_vec=gsl_vector_view_array(pparams+2*S*S+4*S,S);		// body masses
	gsl_vector *Mvec=&M_vec.vector;	

	gsl_matrix_view C_mat=gsl_matrix_view_array(pparams+2*S*S+5*S,S_b,N);		// nutrient contents in plants
	gsl_matrix *Cmat=&C_mat.matrix;											// (P rows, N columns)

	gsl_matrix_view K_mat=gsl_matrix_view_array(pparams+2*S*S+5*S+S_b*N,S_b,N);	// nutrient uptake half saturation densities
	gsl_matrix *Kmat=&K_mat.matrix;											// (P rows, N columns)

	gsl_vector_view S_vec=gsl_vector_view_array(pparams+2*S*S+5*S+2*S_b*N,N);	// nutrient supply concentrations
	gsl_vector *Svec=&S_vec.vector;	

	gsl_vector_view U_vec=gsl_vector_view_array(pparams+2*S*S+5*S+2*S_b*N+N,S_b);// max. nutrient uptake rates
	gsl_vector *Uvec=&U_vec.vector;	

	gsl_vector *Prey=gsl_vector_calloc(S);									// prey biomasses
	gsl_vector *Nvec=gsl_vector_calloc(S);									// net biomass production
	gsl_vector *Dvec=gsl_vector_calloc(S);									// death (mortality+predation losses)

	gsl_matrix *mmat=gsl_matrix_calloc(S,S);
	gsl_vector *rvec=gsl_vector_calloc(S);
	gsl_vector *svec=gsl_vector_calloc(S);
	gsl_vector *nvec=gsl_vector_calloc(N);
	gsl_vector *mvec=gsl_vector_calloc(N);
	gsl_vector *pvec=gsl_vector_calloc(S_b);

//	************ in case of Type 3 functional response ***********
	gsl_vector_memcpy(Prey,Bvec);

	for(i=0;i<S;i++)
	{
		if(gsl_vector_get(Prey,i)<0)
			gsl_vector_set(Prey,i,0);				// paranoia (this doesn't mean that you can safely remove these lines!)

		gsl_vector_set(Prey,i,pow(gsl_vector_get(Prey,i),gsl_vector_get(Hvec,i)));
	}

//	************* functional responses zusammenbauen *************
	gsl_matrix_memcpy(mmat,Hmat); 							//mmat = hij
	gsl_matrix_mul_elements(mmat,Amat);						//mmat = hij*aij

	gsl_blas_dgemv(CblasNoTrans,1,mmat,Prey,0,rvec);		// handling time term: r_i=sum_j(a_ij*h_ij*B_j) 
	gsl_vector_memcpy(svec,Bvec);							// s_i = B_i
	gsl_vector_mul(svec,Cvec);								// predator interference: s_i=C_i*B_i
	gsl_vector_add(rvec,svec);								// r_i = C_i*B_i+sum_j(a_ij*h_ij*B_j)
	gsl_vector_add_constant(rvec,1);						// 1+C_i*B_i+sum_j(a_ij*h_ij*B_j)
	gsl_vector_mul(rvec,Mvec);								// rvec: denominator of functional response r_i = (1+C_i*B_i+sum_j(a_ij*h_ij*B_j) )*M_i

	gsl_blas_dgemv(CblasNoTrans,lambda,Amat,Prey,0,Nvec);	// numerator of prey intake term: N_i=sum_j(A_ij*B_j)

	gsl_vector_div(Nvec,rvec);								// ... dividiert durch den Nenner... N_i = sum_j(A_ij*B_j)/(1+C_i*B_i+sum_j(a_ij*h_ij*B_j) )*M_i
	gsl_vector_sub(Nvec,Rvec);								// ... Respirationsrate abgezogen = (sum_j(A_ij*B_j))/((1+C_i*B_i+sum_j(a_ij*h_ij*B_j))*M_i)-R
	gsl_vector_mul(Nvec,Bvec);								// multiply with target species' biomass  (sum_j(A_ij*B_j)*B_i)/((1+C_i*B_i+sum_j(a_ij*h_ij*B_j))*M_i)  - R*B_i

	gsl_vector_memcpy(svec,Bvec);							// predator biomass  s_i = B_i
	gsl_vector_div(svec,rvec);								// ... mit Nenner der functional response multipliziert s_i = B_i/(1+C_i*B_i+sum_j(a_ij*h_ij*B_j) )*M_i
	gsl_blas_dgemv(CblasTrans,1,Amat,svec,0,Dvec);			// Praedationsverluste pro Individuum: D_i = sum(A_ij * (B_i/(1+C_i*B_i+sum_j(a_ij*h_ij*B_j) )*M_i) )
	gsl_vector_mul(Dvec,Prey);								// multiply with target species' biomass  D_i = sum(A_ij * (B_i/(1+C_i*B_i+sum_j(a_ij*h_ij*B_j) )*M_i) )*B_j

//	************* Terme fuer Basalarten zufuegen *************
	gsl_vector_view u_vec=gsl_vector_subvector(rvec,0,S_b);		// basal species growth rates
	gsl_vector *uvec=&u_vec.vector;
	for(i=0;i<S_b;i++)
	{
		gsl_vector_view K_vec=gsl_matrix_row(Kmat,i);
		gsl_vector *Kvec=&K_vec.vector;							// Kvec is a vector of length N (Kvec_j=Kmat_ij)

		gsl_vector_memcpy(nvec,NBvec);
		gsl_vector_add(nvec,Kvec);
		gsl_vector_memcpy(mvec,NBvec);
		gsl_vector_div(mvec,nvec);								// m_j=NB_j/(Kvec_j+NB_j)

		gsl_vector_set(uvec,i,gsl_vector_min(mvec));			// remember: uvec is just a vector_view of rvec!
	}
	gsl_vector_mul(uvec,Uvec);
	gsl_vector_mul(uvec,PBvec);									// uvec_i=PB_i*mas_i^-0.25*Min_j(NB_j/(NB_j+Kmat_ij))

	gsl_vector_memcpy(nvec,Svec);								// five lines nutrient dynamics...
	gsl_vector_sub(nvec,NBvec);
	gsl_vector_scale(nvec,D_N);
	gsl_vector_memcpy(NBdotvec,nvec);							// NBdot_i=D_N*(S_N_i-N_i)
	gsl_blas_dgemv(CblasTrans,-1,Cmat,uvec,1,NBdotvec);			// NBdot_i=D_N*(S_N_i-N_i)-\sum_j (C_ji*G_j*P_j)

	gsl_vector_mul(rvec,Pvec);									// growth function only for plant species
	gsl_vector_add(Nvec,rvec);

//	************* dynamische Gleichung (Bdot) zusammenfuegen *************
	gsl_vector_memcpy(Bdotvec,Nvec);
	gsl_vector_sub(Bdotvec,Dvec);

	//D_N = 0.25
	//S_vec_j = nutrient supply concentration
	
	sum_j(-t(Cmat_ij*uvec))+dN/dt = D_N(Svec-NBvec)
	
	dN/dt = Sum_j(-t(Cmat_ij)*uvec_j+D_N*(Svec_j-NBvec_j))
	dB/dt = (Nvec+rvec*Pvec)-Dvec
	
//	********** free memory ************
	gsl_vector_free(TBvec);
	gsl_vector_free(Prey);
	gsl_vector_free(Nvec);
	gsl_vector_free(Dvec);
	gsl_vector_free(svec);
	gsl_vector_free(rvec);
	gsl_vector_free(nvec);
	gsl_vector_free(mvec);
	gsl_vector_free(pvec);
	gsl_matrix_free(mmat);

	return GSL_SUCCESS;
}









