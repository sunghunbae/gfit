#include "gfit.h"

using namespace std;

void initialize (CSP &A,int NR,int NX)
{
int i;

/* experimental data */
A.NR = NR;
A.NX = NX;
A.NP = NUMPAR;
A.flag 	= new bool [NR];
A.resid = new unsigned int [NR];
for (i=0;i<NR;i++) {
	A.flag[i] = false;
	A.resid[i] = 0;
	}

A.x   	= gsl_matrix_calloc (NR,NX);
A.y   	= gsl_matrix_calloc (NR,NX);
A.dy 	= gsl_matrix_calloc (NR,NX);

/* parameter */
A.pstat = new unsigned char [NUMPAR]; 
A.p   	= new double [NUMPAR];
A.lb  	= new double [NUMPAR];
A.ub  	= new double [NUMPAR];
A.step 	= new double [NUMPAR];
A.conv 	= new double [NUMPAR];
A.Kl  	= new double [NUMPAR];
A.Ku  	= new double [NUMPAR];
A.mcave = new double [NUMPAR];
A.mcstd = new double [NUMPAR];

for (i=0;i<NUMPAR;i++) {
	RESET_P(A.pstat[i]);
	A.p[i] = A.lb[i] = A.ub[i] = A.step[i] = A.conv[i] = 0.0;
	A.Kl[i] = A.Ku[i] = 1000.0;
	A.mcave[i] = A.mcstd[i] = 0.0;
	}

strcpy(A.pid[_const_],"const");
strcpy(A.pid[_label_],"label");
strcpy(A.pid[_Kd_],"Kd");

/* set defaults */
A.MC = 300;
A.MC_trim = 0;
}

void cleanup (CSP &A)
{
/* free memory */	
delete [] (A.flag);
delete [] (A.resid);
gsl_matrix_free(A.x);
gsl_matrix_free(A.y);
gsl_matrix_free(A.dy);
delete [] (A.pstat);
delete [] (A.p);
delete [] (A.lb);
delete [] (A.ub);
delete [] (A.step);
delete [] (A.conv);
delete [] (A.Kl);
delete [] (A.Ku);
delete [] (A.mcave);
delete [] (A.mcstd);
}

bool converged (CSP &A, double *prev)
/* test if parameters converged */
{
for (int i=0; i< A.NP; i++) {
	if (!IS_FIXED(A.pstat[i]) && fabs(A.p[i]-prev[i]) > A.conv[i]) 
		return false;
	}
return true;
}
