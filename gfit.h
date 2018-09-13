#ifndef _GFIT_H_
#define _GFIT_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <complex>
#include <list>
#include <map>
#include <set>
#include <vector>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>

using namespace std;

#define VERSION	"1.0 (November 2007)"
#define INDENT	4
#define PROMPT	"gf> "

#define NUMPAR	3	// number of parameters
#define LENPAR	20	// length of parameter name

#define LOCAL_KD	0
#define GLOBAL_KD	1

#define _const_		0
#define _label_		1
#define _Kd_            2

/* definitions of pstat */
#define P_RESET	0x00  // 0000 0000
#define P_FIXED	0x01  // 0000 0001 flag: fix(1) or fit(0)
#define P_GLOBL 0x02  // 0000 0010 range: global(1) or local(0)
#define P_UPPER	0x10  // 0001 0000 restrained by upper bound
#define P_LOWER	0x20  // 0010 0000 restrained by lower bound

/* manipulate pstat */
#define RESET_P(x)	(x = P_RESET)
#define SET_FIXED(x)	(x |= P_FIXED)
#define SET_GLOBL(x)	(x |= P_GLOBL)
#define SET_UPPER(x)	(x |= P_UPPER)
#define SET_LOWER(x)	(x |= P_LOWER)
#define UNSET_FIXED(x)  (x &= (~P_FIXED))
#define UNSET_GLOBL(x)  (x &= (~P_GLOBL))
#define UNSET_UPPER(x)  (x &= (~P_UPPER))
#define UNSET_LOWER(x)  (x &= (~P_LOWER))
#define TOGGLE_FIXED(x)	(x ^= P_FIXED)
#define TOGGLE_GLOBL(x)	(x ^= P_GLOBL)
#define TOGGLE_UPPER(x)	(x ^= P_UPPER)
#define TOGGLE_LOWER(x)	(x ^= P_LOWER)

/* evaluate pstat. non-zero if the flag is set */
#define IS_FIXED(x)	(x & P_FIXED)
#define IS_GLOBL(x)	(x & P_GLOBL)
#define IS_UPPER(x)	(x & P_UPPER)
#define IS_LOWER(x)	(x & P_LOWER)

inline double SQR (double a) {return a*a;}

/* chemical shift perturbation */
typedef struct {
	/* data */
	int NR,NX;
	bool *flag;
	unsigned int *resid;
	gsl_matrix *x,*y,*dy;

	/* parameters */
	int NP;
	char pid[NUMPAR][LENPAR];
	unsigned char *pstat;
	double *p,*lb,*Kl,*ub,*Ku,*step,*conv;
	double *mcave, *mcstd;

	/* for minimization */
	unsigned int r,func,mpar,MC, MC_trim;
	double critx2;
	} CSP;

/* func.c */
void 	Kd (void *, gsl_vector *, gsl_matrix *);
double 	chi2 (void *);

/* fitcurve.c */
void 	fitcurve (CSP &);

/* levmar.c */
void 	levmar (CSP &, double &);
int 	chi_f (const gsl_vector *,void *,gsl_vector *);
int 	chi_df (const gsl_vector *,void *,gsl_matrix *);
int 	chi_fdf (const gsl_vector *,void *,gsl_vector *, gsl_matrix *);

/* brent.c */
void 	brent (CSP &, double &);
double 	brent_f (double, void *);

/* gridsearch.c */
double 	grid_search (CSP &);
void 	recursive_grid_search (CSP &, int);
void 	report_grid_search (double, CSP &);

/* file.c */
void 	chk_data (const char *,int &,int &);
void 	read_data (const char *,CSP &);
void 	read_parameter (const char *,CSP &);
void 	parse (const char *,vector <string> &);

/* data.c */
void 	initialize (CSP &,int,int);
void 	cleanup (CSP &);
bool	converged (CSP &, double *);

/* print.c */
void 	print_help ();
void 	terminate (const string);
void 	print_verb (const string);
void 	print_selected (CSP &);
void 	print_iteration (int,CSP &);
void 	print_grid_search (CSP &);
void 	print_boundary (CSP &);

#endif /* _GFIT_H_ */
