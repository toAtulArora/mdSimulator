Molecular Dynamics Simulator
==

The objective of this code will be to understand molecular dynamics by simulating it and verifying it against the results we know analytically. Finally, we might even consider trying non-equilibrium initial points.

[http://www.futa.edu.ng/materials/fortran-novella-Holcomb.pdf]

List of Pre-requisite Objectives
--
	* Make and save animated graphs
	* [todo: complete this]

List of applications you'd need
--
	* synaptic
	* git
	* gfortran
	* gdb
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
	March 27, 2015: Fixed the equation of state! We weren't dividing by 6 in the pressure calculation and we got confused with avgE to mean per particle E, whereas it meant avg over time. We haven't yet figured why in the previous attempt we were getting the same results regardless of collision testing.
	March 30/31, April 1, 2015: Reading about types, classes, elemental functions etc in fortran
	April 2, 2015: Started modernizing. Linked list done, gdb started use
	April 3, 2015: Working on nearest neighbour for a given radius, LJ potential, r-4 potential, gdb used in a big way for debugging, euler evolution, rk4 evolution (neighter worked well with LJ potential, r-4 is fine)
	April 4, 2015: Tweaked rendering for constant physical time, Implemented neighbour search in the real program, further optimization left
	April 5, 2015: (couldn't sleep) (in 2 hours) fixed RK4. Now LJ works fine. able to go upto 1000 particles.
	
Immediate Goals
--
	* Start from non-equilibrium
	* Linked list (perhaps later) | done
	* Neighbouring Boxes | done
	* Short range forces | done (not efficient enough | RK4 was wrong, now fixed, still efficiency can be improved)