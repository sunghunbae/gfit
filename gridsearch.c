#include "gfit.h"

using namespace std;

// definition of grid points
// lb ... ub, step
// 1st grid point = lb + step
// last grid point = ub - step
// (grds - 1) points
// lb and ub are excluded from grid evaluation

double grid_search (CSP &A)  
{
extern int g_neval,g_last;
extern double g_x2min,g_ini[];

g_neval=g_last=0;
g_x2min=-1.0;
for (int i=0;i<A.NP;i++) g_ini[i]=0.0;
recursive_grid_search (A,0);
for (int i=0;i<A.NP;i++) A.p[i] = g_ini[i];
return g_x2min;
}

void report_grid_search (double chisq, CSP &A)
{
extern int g_neval;
printf("%6d ",g_neval);
for (int i=0;i<A.NP;i++)
	if (!IS_FIXED(A.pstat[i])) 
		printf("%s %.4f ",A.pid[i],A.p[i]);
printf(" : %g\n",chisq);
}

void recursive_grid_search (CSP &A,int depth)
{
extern int g_neval,g_last;
extern double g_x2min, g_ini[];

if (depth == 0)
	for(int c=0;c<A.NP;c++) 
		if (!IS_FIXED(A.pstat[c])) g_last=c;

if (depth == A.NP) return;
else { // depth = 0 ... (NP-1)
	if (!IS_FIXED(A.pstat[depth])) {
		A.p[depth] = A.lb[depth];
		while (A.p[depth] <= A.ub[depth]) {
			A.p[depth] += A.step[depth];
			recursive_grid_search (A,depth+1);
			if (depth == g_last) {
				double chisq = chi2 (&A);
				//report_grid_search (chisq, A);
				if (g_x2min == -1.0 || g_x2min > chisq) {
					g_x2min = chisq;
					for (int k=0;k<A.NP;k++) 
						g_ini[k] = A.p[k];
					}
				g_neval++;
				}
			}
		return;
		}
	else recursive_grid_search (A,depth+1);
	};
}
