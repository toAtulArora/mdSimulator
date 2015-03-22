Molecular Dynamics Simulator
==

The objective of this code will be to understand molecular dynamics by simulating it and verifying it against the results we know analytically. Finally, we might even consider trying non-equilibrium initial points.

List of Pre-requisite Objectives
--
	* Make and save animated graphs
	* [todo: complete this]

List of applications you'd need
--
	* synaptic
	* git
	* gfortran
	* gnuplot
	* emacs
	* libavconv

Some Useful Commands
--
	* If you're running ubuntu, might want to do the following if you're running out of space	
	rm -rf ~/.thumbs/*

Progress
--
	March 21, 2015: Attempted to fix equation of state in stage 2. The only possible mistake seems to be that in the initial condition. Attempting to use stage 3 (with collisions turned on) to test. In stage 3 also, we're getting the same result regardless of collision tests. Infact the result doesn't change at all, its same to the last significant digit.