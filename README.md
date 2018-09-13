# gfit
Global Kd Fit

cb 1.0

NAME
	cranberry or shortly cb

SYNOPSIS
	cb [OPTION] [FILE]

DESCRIPTION
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
		assume cb.par is if parameter file is not entered

	-d [DATA FILE]	
		data file. mandatory argument
		see DATA FILE for details

	e.g.
	cb -g -d test.dat or cb -gd test.dat
	cb -l -d test.dat or cb -ld test.dat
	cb -p a.par -gd test.dat
	
DATA FILE 
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
	

PARAMETER FILE
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
	
AUTHOR
	written by Sung-Hun Bae

REPORTING BUGS
	Report bugs to sunghun.bae@gmail.com

COPYRIGHT
	This is free software. You may redistribute copies of it under the
	terms of the GNU General Public License 
	<http://www.gnu.org/licenses/gpl.html> 
	There is NO WARRANTY, to the extent permitted by law.

cb 1.0				November 2007
