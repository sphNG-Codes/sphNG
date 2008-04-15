#!/bin/tcsh
#
# @(#) writes a standard pbs script
#
set phantomdir=~/phantom
echo '#\!/bin/csh'
echo '#PBS -l nodes=1:ppn=8'
echo '#PBS -l walltime=7000:00:00'
echo '#PBS -N '`$phantomdir/scripts/randomword.pl`
echo '#PBS -o pbs.output'
echo '#PBS -j oe'
echo
echo ' cd $PBS_O_WORKDIR'
echo ' echo " PBS_O_WORKDIR $PBS_O_WORKDIR"'
echo ' cat $PBS_NODEFILE >> nodefile'
echo
echo 'echo "starting OpenMP sphNG run"'
echo 'limit stacksize unlimited'
echo 'setenv OMP_SCHEDULE "dynamic,10"'
echo 'setenv OMP_NUM_THREADS 8'
echo
echo './sph_tree_rk_gradh evolution '$1' > '$1'.output'
