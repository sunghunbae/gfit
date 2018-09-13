#include "gfit.h"

using namespace std;

void Kd (void *data, gsl_vector *chi, gsl_matrix *jacob)
{
CSP *A = (CSP *)data;
double t,xi,yi,yexpi,dyi;
double P1 = A->p[_const_];
double P2 = A->p[_label_];
double P3 = A->p[_Kd_];
double fv=0.0;
int i,j;

// impose penalty if boundaries are restrained
for (i=0; i<A->NP; i++) {
	if (!IS_FIXED(A->pstat[i]) && IS_LOWER(A->pstat[i])) 
		fv += A->Kl[i]*SQR(A->p[i] - A->lb[i]);
	if (!IS_FIXED(A->pstat[i]) && IS_UPPER(A->pstat[i])) 
		fv += A->Ku[i]*SQR(A->p[i] - A->ub[i]);
	}

for (i=0; i<A->NX; i++) {
	xi = gsl_matrix_get (A->x,A->r,i);
	t = sqrt(SQR(xi+P2+P3)-(4*xi*P2));
        dyi = gsl_matrix_get (A->dy,A->r,i);
	if (chi != NULL) {
                yi = P1*((xi+P2+P3)-t);
                yexpi = gsl_matrix_get (A->y,A->r,i);
		gsl_vector_set(chi, i, fv + (yi-yexpi)/dyi);
		}
	if (jacob != NULL) {
		j=0;
		if(!IS_FIXED(A->pstat[_const_]))
			gsl_matrix_set(jacob,i,j++,((xi+P2+P3)-t)/dyi);
		if(!IS_FIXED(A->pstat[_label_]))
			gsl_matrix_set(jacob,i,j++,P1*(1-1/t*(P2+P3-xi))/dyi);
		if(!IS_FIXED(A->pstat[_Kd_]))
			gsl_matrix_set(jacob,i,j++,P1*(1-1/t*(P2+P3+xi))/dyi);
		}
	}
}

double chi2 (void *data)
{
CSP *A = (CSP *) data;
double chisq,chisq_grid,X2;
int i,r;

gsl_vector * chi = gsl_vector_calloc (A->NX);

if (A->func == LOCAL_KD) {
	Kd (data, chi, NULL);
	for (X2=0.0,i=0;i<A->NX;i++)
		X2 += SQR(gsl_vector_get(chi,i));
	}

if (A->func == GLOBAL_KD) {

	/* hold global parameter */
	for (i=0;i<A->NP;i++)
		if (IS_GLOBL(A->pstat[i])) TOGGLE_FIXED(A->pstat[i]);

	/* local fitting while holding global parameter */
	A->func = LOCAL_KD;
	A->MC = 0;
	for (X2=0.0,r=0;r<A->NR;r++) {
		A->r = r;
		chisq_grid = grid_search(*A);
		levmar(*A,chisq);
		X2 += chisq;
		}

	/* release global parameter */
	for (i=0;i<A->NP;i++)
		if (IS_GLOBL(A->pstat[i])) TOGGLE_FIXED(A->pstat[i]);

	/* restore function */
	A->func = GLOBAL_KD;
	}

gsl_vector_free (chi);
return X2;
}
