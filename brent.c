#include "gfit.h"
#include <gsl/gsl_min.h>

#define BRENT_ITER_MAX		1000000

using namespace std;

/* 1D global fitting */
/* interval [lower...upper], initial guess */
void brent (CSP &A, double &chisq)
{
double lower,upper,guess,conv,fret;
int i,dim=0,px,status,iter = 0;

for (i=0;i<A.NP;i++)
	if (!IS_FIXED(A.pstat[i]) && IS_GLOBL(A.pstat[i])) {
		dim++;
		px=i;
		}

if (dim > 1) 
	printf("Brent Called for Multi-dimensional Minimization\n");

A.mpar = px;
	
lower=A.lb[px];
upper=A.ub[px];
guess=A.p[px];
conv =A.conv[px];

const gsl_min_fminimizer_type *T = gsl_min_fminimizer_brent;
gsl_min_fminimizer *s = gsl_min_fminimizer_alloc (T);
gsl_function F = {&brent_f, &A};
gsl_min_fminimizer_set (s, &F, guess, lower, upper);

do {
	iter++;
	status = gsl_min_fminimizer_iterate (s);
	fret = gsl_min_fminimizer_f_minimum (s);
	guess = gsl_min_fminimizer_x_minimum (s);
	lower = gsl_min_fminimizer_x_lower (s);
	upper = gsl_min_fminimizer_x_upper (s);
	status = gsl_min_test_interval (lower,upper,conv,0.0);
	if (status == GSL_SUCCESS) A.p[px] = guess;
	printf ("%5d [%.7f, %.7f] %.7f %.7f f= %g\n",
		iter,lower,upper,guess,upper-lower,fret);
	}
while (status == GSL_CONTINUE && iter < BRENT_ITER_MAX);

gsl_min_fminimizer_free (s);
chisq = fret;
}

double brent_f (double x, void * params)
{
CSP *A = (CSP *)params;
A->p[A->mpar] = x;
return chi2 (A);
}
