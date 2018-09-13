#include "gfit.h"
#include <getopt.h>

#define MAX_ITER 100

using namespace std;

/* internal use: seed for random number generation */
long unsigned int seed = time(0);

/* for recursive grid search routine */
int 	g_neval,g_last;
double 	g_ini[NUMPAR];
double 	g_x2min,g_PROB;

int main(int argc, char *argv[])
{
/* signature */
printf("\n\tgfit ");printf(VERSION); 
printf(" written by Sung-Hun Bae\n\n"); 

/* checking elapsed time */
time_t start,finish;
time(&start);

int i,j;

char *D_FILE = NULL;
char *P_FILE = NULL;

/* COMMANDS and OPTIONS */
bool cmd_help = false;
bool cmd_global = false;
bool cmd_local = false;
bool cmd_test = false; // reserved for test 

/* commands arguments */
int c;
opterr = 0;
while ((c = getopt (argc, argv, "hvgltp:d:")) != -1)
	switch (c) {
		case 'h': cmd_help = true; break;
		case 'g': cmd_global = true; break;
		case 'l': cmd_local = true; break;
		case 't': cmd_test = true; break;
		case 'p': P_FILE = optarg; break;
		case 'd': D_FILE = optarg; break;
		case '?': cmd_help = true; break;
           	}// switch

if (cmd_help) print_help();

/* declaration of variables */
const int NP = NUMPAR;
int NR;// number of residues
int NX;// number of data points
int r;

/* check data file & Allocate memory */
print_verb ("check data file");
chk_data (D_FILE,NR,NX); 

/* INITIALIZE MEMORY */
CSP A;
initialize (A,NR,NX);

/* READ PARAMETER FILE */
printf(PROMPT);
printf("reading parameter file : ");
if (P_FILE == NULL) {
	printf("gfit.par\n");
	read_parameter("./gfit.par",A);
	}
else {
	printf("%s\n",P_FILE);
	read_parameter(P_FILE,A);
	}
print_boundary(A);


/* READ DATA FILE */
printf(PROMPT);
printf("reading data file : %s\n",D_FILE);
read_data(D_FILE,A);

/* FIT GLOBAL Kd */
if (cmd_global) {
	double chisq,max_Kd=-1.0,min_Kd=-1.0;
	double X2,chisq_grid;

	print_verb("fit GLOBAL Kd");

	print_verb("first, fitting local Kd");
	for (r=0;r<A.NR;r++) {
		A.r = r;
		A.func = LOCAL_KD;
		A.MC = 300;
		fitcurve (A);
		string indent(INDENT,' ');
		printf("%s",indent.c_str());
		printf("RESIDUE %4d %5s %8.4g +/- %8.4g %5s %8.2e +/- %8.2e\n",
        		A.resid[r],A.pid[_const_],A.p[_const_],A.mcstd[_const_],
        		A.pid[_Kd_],A.p[_Kd_],A.mcstd[_Kd_]);

		if (A.p[_Kd_] < min_Kd || min_Kd == -1.0 ) min_Kd = A.p[_Kd_];
		if (A.p[_Kd_] > max_Kd || max_Kd == -1.0 ) max_Kd = A.p[_Kd_];

		}// r

	SET_GLOBL(A.pstat[_Kd_]);
	A.func = GLOBAL_KD;
	A.p[_Kd_] = 0.5*(min_Kd + max_Kd);
	A.lb[_Kd_] = min_Kd*0.5;
	A.ub[_Kd_] = max_Kd*2.0;
	A.step[_Kd_] = (max_Kd - min_Kd)/50;

	print_verb("set parameter boundary");
	print_boundary(A);

	print_verb("global minimization by BRENT algorithm");
	brent (A, chisq);

	print_verb("adjusted local parameters");

        /* hold global parameter */
        for (i=0;i<A.NP;i++)
                if (IS_GLOBL(A.pstat[i])) SET_FIXED(A.pstat[i]);

        /* local fitting while holding global parameter */
        A.func = LOCAL_KD;
        A.MC = 300;
        for (X2=0.0,r=0;r<A.NR;r++) {
                A.r = r;
		fitcurve(A);
		string indent(INDENT,' ');
		printf("%s",indent.c_str());
                printf("RESIDUE %4d %5s %8.4g +/- %8.4g %5s %8.2e\n",
                        A.resid[r],A.pid[_const_],A.p[_const_],
                        A.mcstd[_const_],A.pid[_Kd_],A.p[_Kd_]);
                }

        /* release global parameter */
        for (i=0;i<A.NP;i++)
                if (IS_GLOBL(A.pstat[i])) UNSET_FIXED(A.pstat[i]);
	A.func = GLOBAL_KD;

	}// cmd_global_model


/* FIT LOCAL Kd */
if (cmd_local) {
	print_verb("FIT LOCAL Kd with Monte Carlo simulations");
	for (r=0;r<A.NR;r++) {
		A.r = r;
		A.func = LOCAL_KD;
		fitcurve (A);
		string indent(INDENT,' ');
		printf("%s",indent.c_str());
		printf("RESIDUE %4d %5s %8.4g +/- %8.4g %5s %8.2e +/- %8.2e\n",
        		A.resid[r],A.pid[_const_],A.p[_const_],A.mcstd[_const_],
        		A.pid[_Kd_],A.p[_Kd_],A.mcstd[_Kd_]);
		}// r
	}// cmd_local_model


/* RESERVED FOR TEST */
if (cmd_test) { }

/* EXIT */
cleanup (A); // free memory

time(&finish);
int hh=0,mm=0,ss=0,elapsed=finish-start;
hh = (int)floor((double)(elapsed/3600));
mm = (int)floor((double)((elapsed/60)-(hh*60)));
ss = (elapsed - (hh*3600)-(mm*60));
printf(PROMPT); printf("started  %s",ctime(&start));
printf(PROMPT); printf("finished %s",ctime(&finish));
printf(PROMPT); printf("time elapsed %02d:%02d:%02d\n",hh,mm,ss);
}
