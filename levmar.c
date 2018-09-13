#include "gfit.h"
#include <gsl/gsl_multifit_nlin.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_version.h>

#define LEVMAR_ITER_MAX         1000000
#define LEVMAR_STOP_GRAD        1E-15
#define LEVMAR_STOP_DELTA       1E-09

/* function for fitting */
/* funciton should provide chi(i) vector and jacob(i,j) matrix */
/* i: number of x or y data points */
/* j: number of active fitting paramters */

void (*func) (void *data, gsl_vector *chi, gsl_matrix *jacob) = &Kd;


using namespace std;

/* Levenberg-Marquardt Nonlinear Least Square Fitting */
void levmar (CSP &A, double &chisq)
{
int i,j,mfit;

for (mfit=0,i=0;i<A.NP;i++) 
	if (!IS_FIXED(A.pstat[i])) mfit++;

gsl_vector * p = gsl_vector_calloc(mfit);
gsl_matrix * covar = gsl_matrix_alloc (mfit,mfit);
gsl_multifit_function_fdf f;
gsl_multifit_fdfsolver *s;
const gsl_multifit_fdfsolver_type *T;
int status;
unsigned int iter = 0;

for (i=0,j=0;i<A.NP;i++)
	if (!IS_FIXED(A.pstat[i])) 
		gsl_vector_set(p,j++,A.p[i]);

f.f = &chi_f;
f.df = &chi_df;
f.fdf = &chi_fdf;
f.n = A.NX;
f.p = mfit;
f.params = &A;
	
T = gsl_multifit_fdfsolver_lmsder;
s = gsl_multifit_fdfsolver_alloc (T, A.NX, mfit);
gsl_multifit_fdfsolver_set (s, &f, p);

// fix GSL2.x compiling error
#if GSL_MAJOR_VERSION > 1
gsl_matrix *J = gsl_matrix_alloc(f.n, f.p);
#endif

do {
	iter++;
	status = gsl_multifit_fdfsolver_iterate (s);
	if (status) break;
	status = gsl_multifit_test_delta (s->dx, s->x, 
		LEVMAR_STOP_GRAD, LEVMAR_STOP_DELTA);
	} 
while (status == GSL_CONTINUE && iter < LEVMAR_ITER_MAX);

chisq = gsl_blas_dnrm2(s->f);
chisq = SQR(chisq);

// fix GSL2.x compiling error
#if GSL_MAJOR_VERSION > 1
gsl_multifit_fdfsolver_jac(s,J);
gsl_multifit_covar(J, 0.0, covar);
#else
gsl_multifit_covar(s->J, 0.0, covar);
#endif

gsl_multifit_fdfsolver_free(s);
gsl_matrix_free (covar);
gsl_vector_free (p);
}


int chi_f (const gsl_vector *p, void *data, gsl_vector *chi)
/* function definition : f(i) = (y(i)-yexp(i))/dy(i) */
/* y(i): calculated value for ith data */
/* yexp(i): experimental value for ith data */
/* dy(i): uncertainty for ith data */
{
CSP *A = (CSP *)data;
int i,j;

// project parameters
for (i=0,j=0; i<A->NP; i++)
	if (!IS_FIXED(A->pstat[i]))
		A->p[i] = gsl_vector_get (p,j++);

func (data,chi,NULL);

return GSL_SUCCESS;
}

int chi_df (const gsl_vector *p, void *data, gsl_matrix *jacob)
/* jacobian definition : J(ij) = df(i)/dp(j) */
/* df(i)/dp(j): derivative of f(i) with respect to p(j) parameter */
{
CSP *A = (CSP *)data;
int i,j;

// project parameters
for (i=0,j=0; i<A->NP; i++)
	if (!IS_FIXED(A->pstat[i]))
		A->p[i] = gsl_vector_get (p,j++);

func (data,NULL,jacob);
return GSL_SUCCESS;
}

int chi_fdf (const gsl_vector *p, void *data, gsl_vector *chi, 
	gsl_matrix *jacob)
/* returns function value and jacobian matrix */
{
CSP *A = (CSP *)data;
int i,j;

// project parameters
for (i=0,j=0; i<A->NP; i++)
	if (!IS_FIXED(A->pstat[i]))
		A->p[i] = gsl_vector_get (p,j++);

func (data,chi,jacob);
return GSL_SUCCESS;
}
