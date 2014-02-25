This repository contains MATLAB and Python scripts for metabolic network analysis and  performing Flux Balance Analysis (FBA). The repository is hosted privately at BitBucket and mirrored publicly at GitHub.


File and resource types:
------------------------

**.m files:** MATLAB script files. May be opened/executed with MATLAB or GNU Octave.

**.mat files:** MATLAB data files. To open a .mat file in GNU Octave, use the load command:
"load <filename>"

**.py files:** Python  2.7 script files.

MATLAB files:
-------------

The MATLAB scripts were developed for a half-semester project at the Norwegian University of Scienc and Technology (NTNU) in fall 2012. The scripts can be used to compare the results of Flux Balance Analysis (FBA) with experimentally determined flux data, calculate the euclidean distance between experimental and computed fluxes and determine which model parameters (combinations of objective function and constraints) which are most consistent with the experimental data. The MATLAB scripts will not be developed further, as they are superceded by the Python scripts.

Python files:
-------------

The collection of Python scripts expand on the functionality of the MATLAB files, and also include implementations of several non-related methods previously only available as MATLAB scripts or not available as open source code. The scripts were written starting in in Fall 2013. Most scripts requires the free FBA package CobraPy and a linear optimization solver supported by that package, such as Gurobi. NumPy and/or SciPy are also required for most of the scripts.

A description of some of the scripts and functions follows below:

**optreqanalysis.py:**
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

**QPmindist.py:**

Uses quadratic programming in conjunction with flux balance analysis (FBA) to determine the minimal euclidean distance between a set of flux values and the solution of an FBA problem, for a given
requirement of objective-optimality for the FBA solution.

	Dependencies:
		compdist.py (test block only)
		extractflux.py
		Gurobipy
		Cobrapy

	Reads from:
		reactionmaps.mat (test block only)


**Compdist.py:**

Loads flux values from stored experimental data and compares with a supplied flux vector (the first argument). Returns the euclidean distance between the two vectors.

	Dependencies:
	Experimental data.
	loadData.py
	extractflux.py


**constrainfluxes.py:**

The function constrainfluxes takes as input a COBRA model and constrains one or more reactions with a dened corresponding experimental reaction to limits specified by the experimental value and by a tolerance relative to the
experimental uncertainty.

Which fluxes to constrain is specified by the optional field "fluxconstrainset", which is a list containing the reactions to be constrained. By default, all fluxes for which experimental data is available are constrained. Default relative tolerance is 1.

    Input:
    	*A CobraPy model
    	*Experimental fluxes
    	*experimental errors
    	*reaction map
    	*error tolerance (optional)
    	*fluxconstrainset (optional)
    	
    Output:
    	*Modified CobraPy model

	
**extractflux.py:**
The function extractflux takes as input a flux vector from FBA result and returns a subset of those reactions as a list or dictionary with fluxes corresponding to the experimental reactions.

	Used by:
	QPmindist.

**Convertrxmaps:**

Loads legacy reaction maps from a .mat file.

**exportreactionIDs.py:**

Writes all reaction IDs in a CobraPy model to a text file.

**Fluxreport.py:**

The function fluxreport takes as input a vector of fluxes corresponding to the set of experimental reaction rates, and returns both sets of fluxes, the difference between each corresponding flux, the uncertainty in the experimental fluxes and the difference between the fluxes divided by that uncertainty.


Data files:
-----------
For the scripts comparing experimental and computed flux values, the following data must be supplied :

**Experimental data:**
Experimental flux values must be supplied in .xml format. See expdata.xml for example data.

**Reaction maps:**

Reaction maps are loaded from  .txt/ or .xml files by the script loadData.py. See reactionmaps.xml for an example file.





