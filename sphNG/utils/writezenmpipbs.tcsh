#!/bin/tcsh
#
# @(#) writes a standard pbs script
#
if $# != 1 then
   echo Usage: $0 infile
else
   set phantomdir=~/phantom
   echo '#\!/bin/csh'
   echo '#PBS -l nodes=8:ppn=8:StandardMem'
   echo '#PBS -l walltime=7000:00:00'
   echo '#PBS -N '`$phantomdir/scripts/randomword.pl`
   echo '#PBS -o pbs.output'
   echo '#PBS -j oe'
   echo
   echo 'cd $PBS_O_WORKDIR'
   echo 'echo " PBS_O_WORKDIR $PBS_O_WORKDIR"'
   echo 'cat $PBS_NODEFILE'
   echo
   echo 'echo "starting OpenMP+MPI sphNG run"'
   echo 'limit stacksize unlimited'
   echo 'setenv OMP_SCHEDULE "dynamic,10"'
   echo 'setenv OMP_NUM_THREADS 8'
   echo
   echo '# for Intel MPI'
   echo 'setenv I_MPI_DEVICE rdssm'
   echo '#setenv I_MPI_PIN_MODE lib'
   echo '#setenv I_MPI_PIN_PROCS 0,2'
   echo '#setenv I_MPI_DEBUG 3'
   echo
   echo 'setenv NP `cat ${PBS_NODEFILE}|wc -l`'
   echo 'cat ${PBS_NODEFILE} > tmp'
   echo 'sort -u tmp > mpd.hosts'
   echo 'setenv NODES `wc -l < mpd.hosts`'
   echo 'mpdboot -n ${NODES} -r ssh --file=mpd.hosts'
   echo 'echo ${NP}'
   echo 'mpiexec -perhost 1 -envall -np ${NODES} dplace -x3 ./sph_tree_rk_gradh_RT evolution '$1' > '$1'.output'
   echo
   echo 'mpdallexit'
endif
