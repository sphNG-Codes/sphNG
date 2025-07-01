# sphNG

Project Summary:

The main aim of this project is to add GPU support to the astrophysical Smoothed Particle Hydrodynamics (SPH) code SPHNG. The code is written in Fortran and currently parallelised using a hybrid MPI-OpenMP approach.

sphNG https://bitbucket.org/tb-astro/sphng/src/master/

Thomas Bending email: t.j.r.bending@exeter.ac.uk 

# Documentation

# Installation

# Codes

A. GPU implementation of the radiative transfer (RT) calculations. The most expensive OpenMP loop in the code calculates RT using the flux-limited diffusion approximation and using an iterative solver.
B. Neighbour finding.
C. Standalone SPH density solver.

# Problems 


