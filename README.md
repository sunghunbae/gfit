# gfit 1.0 (November 2007)
Global Kd Fit

## SYNOPSIS
	gfit [OPTION] [FILE]

## DESCRIPTION
```
	Performs optimization and non-linear least square fitting of 
	a function that describes the relation between experimental 
	parameters including dissociation constant Kd and measured NMR
	chemical shifts perturbation data. During fitting, Kd can be 
	treated as a global or local parameter. By definition, the value 
	of a global parameter is shared by all selected residues whereas 
	that of local parameter is not shared and only valid for one residue.
	
	-l	(default)		
		local fit. selected when neither -g nor -l is given 

	-g		
		global fit

	-p [PARAMETER FILE]	
		parameter file.
		see PARAMETER FILE for details
		assume gfit.par is if parameter file is not entered

	-d [DATA FILE]	
		data file. mandatory argument
		see DATA FILE for details

	e.g.
	gfit -g -d test.dat or gfit -gd test.dat
	gfit -l -d test.dat or gfit -ld test.dat
	gfit -p a.par -gd test.dat
```
	
## DATA FILE 
```
	Single data file in ASCII text format is required.
	Lines begining with '#' are regared as comments.
	The first line of data set for a residue begins by keyword 'RESIDUE'.
	Residue number and flag('1' or '0') for fitting follow 'RESIDUE'.
	The residue with flag 0 will be excluded in the fitting. 
	From the next line, data lines follow.
	Each data line has X Y DY format in which DY means uncertainty of Y.
	X, Y, and DY should be seperated by at least single space.
	The last line of data set for a residue is remarked by '&'

	e.g.
	RESIDUE 121 1
	0.010	0.023	0.01
	0.025	0.042	0.01
	.....	.....	....
	0.200	0.314	0.01
	0.250	0.330	0.01
	0.300	0.335	0.01
	&
	RESIDUE 122 0
	.....	.....	.....
```
	
## PARAMETER FILE
```
	Single parameter file in ASCII text format is required.
	Reserved keywords for the first words: 'fit', 'par', 'fix'
	Reserved parameter names: 'const', 'label', 'Kd'

	fit -mc [INTEGER] -mctrim [INTEGER] -critx2 [REAL]
	: define Monte Carlo simulations and chi-square statistics for
	 estimation of fitting error and acceptance criteria, respectively. 

	par [PARAMETER_NAME] [LOWER_BOUND] [UPPER_BOUND] [GRID_STEP] 	\
		[CONVERGENCY_LIMIT] [PENALTY_IF_LOWER_BOUND_VIOLATED] 	\
		[PENALTY_IF_UPPER_BOUND_VIOLATED]
	: define boundaries for grid search and restrained fitting

	fix [PARAMETER_NANE] [REAL]
	: fix given parameter to the given value during fitting.
	
	e.g.
	fit -mc 300 -mctrim 5 -critx2 0.05
	par const 0 10 1 0.001 0 0
	par Kd 0 0.1 0.001 0.0001 0 0
	fix label 0.2
```

## Example Output

```

	gfit 1.0 (November 2007) written by Sung-Hun Bae

gf> check data file
    name of data file                   : test.dat
    number of lines                     : 93
    number of residues to be fitted     : 3 out of 3
    number of maximum data points       : 29
gf> reading parameter file : test.par
    number of MC simulations 300
    number of MC simulations trimmed 2*5
    acceptance for Chi-square criteria 0.05
    fix label = 0.2
       Param     Value [ Lower BND ... Upper BND | Grid Step ] Conv.Limit
    L* const 0.000e+00 [ 0.000e+00 ... 1.000e+01 | 1.000e+00 ] 1.00E-03
    L  label 2.000e-01 [ 0.000e+00 ... 1.000e+00 | 1.000e-02 ] 1.00E-02
    L* Kd    0.000e+00 [ 0.000e+00 ... 1.000e-01 | 1.000e-03 ] 1.00E-04
gf> reading data file : test.dat
gf> fit GLOBAL Kd
gf> first, fitting local Kd
    RESIDUE    1 const      2.5 +/-    0.013    Kd 9.99e-03 +/- 8.00e-04
    RESIDUE    2 const        6 +/-  0.01392    Kd 1.50e-02 +/- 4.34e-04
    RESIDUE    3 const     15.5 +/-  0.01409    Kd 1.30e-02 +/- 1.52e-04
gf> set parameter boundary
       Param     Value [ Lower BND ... Upper BND | Grid Step ] Conv.Limit
    L* const 1.550e+01 [ 0.000e+00 ... 1.000e+01 | 1.000e+00 ] 1.00E-03
    L  label 2.000e-01 [ 0.000e+00 ... 1.000e+00 | 1.000e-02 ] 1.00E-02
    G* Kd    1.249e-02 [ 4.994e-03 ... 2.998e-02 | 1.000e-04 ] 1.00E-04
gf> global minimization by BRENT algorithm
    1 [0.0049940, 0.0191681] 0.0124881 0.0141741 f= 54.647
    2 [0.0124881, 0.0191681] 0.0130861 0.0066800 f= 32.708
    3 [0.0130861, 0.0191681] 0.0131831 0.0060820 f= 32.4757
    4 [0.0130861, 0.0131831] 0.0131589 0.0000970 f= 32.4496
gf> adjusted local parameters
    RESIDUE    1 const    2.545 +/- 0.005766    Kd 1.32e-02
    RESIDUE    2 const    5.943 +/- 0.006078    Kd 1.32e-02
    RESIDUE    3 const    15.51 +/- 0.005937    Kd 1.32e-02
gf> started  Thu Sep 13 10:48:11 2018
gf> finished Thu Sep 13 10:48:12 2018
gf> time elapsed 00:00:01
```
