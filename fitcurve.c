#include "gfit.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

using namespace std;

void fitcurve (CSP &A)
{
int i,j,k;
double chisq_grid, chisq;

/* degree of freedom */

/* do grid search for starting point */
chisq_grid = grid_search (A);

/* Levenberg-Marquardt Nonlinear Least Square Fitting */
levmar (A,chisq);

if (chisq >= chisq_grid) {
	printf("WARNING: Nonlinear Least Square fitting may be failed\n");
	}
	
if (A.MC > 0) {
	// number of simulations
	int mc_trial=0, mc_success=0;

	double exp[A.NX],par[A.NP];

	// backup experimental data
	for (i=0;i<A.NX;i++)
		exp[i] = gsl_matrix_get(A.y,A.r,i);

	// backup fitted parameters
	for (i=0;i<A.NP;i++)
		par[i] = A.p[i];

	// random gaussian distribution
	extern long unsigned int seed;
	const gsl_rng_type * T = gsl_rng_default;
	gsl_rng * random = gsl_rng_alloc (T);
	gsl_rng_set (random,seed);
	double sim;

	// memory allocation for simulated data
	gsl_matrix *mcdat = gsl_matrix_calloc(A.MC,A.NP);

	// error estimation
	double s,ave,var,ep;
	int n;

	/* Monte Carlo Simulations for error estimation *************/
	do {
		for(j=0;j<A.NX;j++) {
			sim = gsl_ran_gaussian(random,
				gsl_matrix_get(A.dy,A.r,j));
            		gsl_matrix_set(A.y, A.r, j, exp[j] + sim); 
			}

		/* do grid search for new starting point */
		chisq_grid = grid_search (A);

		/* Levenberg-Marquardt Nonlinear Least Square Fitting */
		levmar (A,chisq);
		if (chisq < chisq_grid) {
			for (i=0;i<A.NP;i++) 
				gsl_matrix_set(mcdat,mc_success,i,A.p[i]);
			mc_success++;
			} // account only successful fitting 

		mc_trial++;

		if (mc_trial > A.MC*4) break;

		} while (mc_success < A.MC); // MC loop 
	
	/* restore experimental data */
	for (i=0;i<A.NX;i++)
		gsl_matrix_set(A.y,A.r,i,exp[i]);

	/* restore experimental data */
	for (i=0;i<A.NP;i++)
		A.p[i] = par[i];
		
	/* Estimate error as standard deviation of fitted values */
	for (i=0;i<A.NP;i++) {
		s=0.0;
        	for (n=0,j=A.MC_trim;j<(mc_success-A.MC_trim);j++) {
            		s += gsl_matrix_get(mcdat,j,i);
			n++;
			}
		A.mcave[i] = ave = s/n;

        	ep=var=0.0;
        	for (j=A.MC_trim;j<(mc_success-A.MC_trim);j++) {
            		s = gsl_matrix_get(mcdat,j,i)-ave;
            		ep += s;
            		var += s*s;
            		}
        	var=(var-ep*ep/n)/n;
		A.mcstd[i] = sqrt(var); // stdev
        	} // i

	gsl_matrix_free (mcdat);

	// change seed number for random number generator
	seed = (long unsigned int)gsl_rng_get(random); 
	gsl_rng_free (random);

	} // if MC > 0
}
