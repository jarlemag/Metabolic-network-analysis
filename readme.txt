This repository contains various MATLAB and Python scripts for analyzing metabolic networks and performing Flux Balance Analysis (FBA). 



MATLAB files:

The MATLAB version was developed for a half-semester project at NTNU fall 2012. The scripts can be used to compare the results of Flux Balance Analysis (FBA) with experimentally determined flux data, calculate the euclidean distance between experimental and computed fluxes and determine which model parameters (combinations of objective function and constraints) which are most consistent with the experimental data. The MATLAB scripts will not be developed further, as they are being superceded by the Python scripts.



Python:

The Python scripts expand on the functionality of the MATLAB files, and also include implementations of several non-related methods previously only available as MATLAB scripts or not available as open source code. The scripts were written starting in in Fall 2013. Most scripts requires the free FBA package CobraPy, a linear optimization solver supported by that package, such as Gurobi. NumPy and/or SciPy are also required for most of the scripts.

A description of some of the scripts and functions follows below:

optreqanalysis.py:
Used to evaluate the minimal achievable distance between a FBA solution and a set of experimental fluxes, for various objective-optimality requirements with respect to the FBA solution.
Can process several models and conditions (objectives, constraints) in one run. 
Receives as first argument a Python dictionary containing one or several CobraPy models, along with a reaction map for each model that specifies which reactions in the model correspond to which reactions in the experimental data. 
The third argument is a list of experiments for which experimental data exists, and for which the computations should be made.
The second argument is a list of objective functions which should be used. 

Loads one or several reaction maps relating reactions and flux values in an experimental data set to appropriate reactions in the model(s) to be used.



	Dependencies:

	QPmindist.py

	
	Reads from:
		expdata.mat
	Writes to:


QPmindist.py:

Uses quadratic programming in conjunction with flux balance analysis (FBA) to determine the minimal euclidean distance between a set of flux values and the solution of an FBA problem, for a given
requirement of objective-optimality for the FBA solution.



	Dependencies:
		compdist2.py (test block only)
		extractflux2.py
		Gurobipy
		Cobrapy

	Reads from:
		reactionmaps.mat (test block only)



Compdist:
Obsolete version of compdist2 based on MATLAB code, but with more options. Should consider adding options from compdist to compdist2.


Compdist2.py:

Loads a vector of flux values from stored experimental data and compares with a supplied flux vector (the first argument).
Returns the euclidean distance between the two vectors.


	Dependencies:
	Experimental data.		

	Reads from:
	expdata.mat (currently), should change this to JSON, XML or .txt/csv.

	Current issues: Experimental data is hard-coded to Perrenoud batch aerobe data set. Need to update script to accept arbitrary experimental data.

constrainfluxes.py:

The function constrainuxes takes as input a COBRA model and constrains
one or more reactions with a dened corresponding experimental reaction to
limits specied by the experimental value and by a tolerance relative to the
experimental uncertainty.

Input: 
	*A CobraPy model
	*An options dictionary


Output:
	*Modified CobraPy model

Which fluxes to constrain is specified by the optional field "fluxconstrainset", which is a list containing the reactions to be constrained.
As default, all fluxes (for which experimental data is available) are constrained. 

Default relative tolerance is 1.


	Dependencies:

		

	Status: Not operational. Need to make code for use of reaction maps and application of constraints.


In the MATLAB version, constrainfluxes.m is called by runsim.m 
	


extractflux.py:
Obsolete version of extractflux2.py
	
	Dependencies:


	Reads from:
	reactionmaps.mat


extractflux2.py:
The function extractux takes as input a flux vector from FBA result and re-
turns a subset of those reactions as a list or dictionary with fluxes corresponding to the experimental reactions.



	Used by:
	QPmindist.
	
	Current issues:
	Update script to work with dictionaries in place of flux vectors.



Convertrxmaps:

Loads legacy reaction maps from a .mat file.



exportreactionIDs.py:

Writes all reaction IDs in a CobraPy model to a text file.




Data files:


The following files must be supplied:


Experimental data:





Reaction maps:

Reaction maps are loaded from .mat (currently) or .txt/.xml/.JSON/.SBML files (planned) in the script optreqanalysis.py








Fluxreport.py:

The function uxreport takes as input a vector of uxes corresponding to the
set of experimental reaction rates, and returns both sets of uxes, the dier-
ence between each corresponding ux, the uncertainty in the experimental
uxes and the dierence between the uxes divided by that uncertainty.




	Status: Operational

	Current issues: Should add readout on solver, model, objective and constraints.

Example use:







Files and resources:

.m files: MATLAB script files. May be opened/executed with MATLAB or GNU Octave.

.mat files: MATLAB data files. To open a .mat file in GNU Octave, use the load command:
"load <filename>"



.py files: Python  3 script files.



