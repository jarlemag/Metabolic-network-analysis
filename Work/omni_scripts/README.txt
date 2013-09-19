Description:

Matlab scripts and GAMS input files demonstrating the application of the OMNI
method.

Usage:

Run omni_example.m in Matlab

Requirements:

Requires that GAMS and the CPLEX MILP solver are installed and accessible
from this directory.  

Contents:

omni_example.m	Sets up an OMNI problem applied to exchange flux data and
		growth rate data for an experimentally evolved E. coli strain
run_omni.m	Performs preprocessing to set up the actual constraints and
		variables for the OMNI problem
omni.m		Creates the matrices, lower, and upper bounds needed by the
		optimization software (in input format suitable for the LINDO
		MILP solver)
create_gams_input.m  Converts the problem into a form that can be solved by
		     the CPLEX solver through GAMS
print_gams_input.m   Prints the GAMS input file
parse_gams_omni_output.m    Parses the GAMS OMNI solution

omni_driver.gms		    GAMS driver to solve general OMNI problems (in
			    fact to solve any MILP problem)
cplex.opt		    Options for the CPLEX solver

lib_files		    A nubmer of general usage Matlab scripts needed by
			    the main OMNI functions

Disclaimer: This is very much a beta-version. It has been primarily tested
	    with the E. coli data described in the PLoS Comp Biol paper.

Markus Herrgard 2/8/06
mherrgar@ucsd.edu
