#!/bin/bash
#
# script to run a sequence of runs via xgrid
#
if [ "$#" -ne 1 ]; then
  echo Submits multiple runs of sphNG to xgrid server
  echo by D. Price Sep 2006
  echo
  echo Usage: $0 jobtype;
  echo "where job type is 'initial', 'restart' or 'dummy'";
else
jobtype=$1;
if [ $MACHTYPE == 'i386' ]; then
   echo 'INTEL MAC GRID'
   gid=1;
else
   echo 'PPC MAC GRID'
   gid=0;
fi
xgridauth='-hostname cytosine.ex.ac.uk -auth Kerberos -gid '$gid;

if [ $jobtype == 'initial' ]; then
   echo initial run;
elif [ $jobtype == 'dummy' ]; then
   echo dummy run;
fi
if [ -f cleanup ]; then
   source cleanup
   rm cleanup
fi
rm getresults
rm checkstatus
##machines=( `cat machinelist` );
for x in `cat masstoflux_values`; do 
newdir=mbossbod_f$x
mkdir $newdir
cd $newdir
echo '--------' $newdir '--------'
pwd
#
# make sure we have the latest code version
#
~/sphNG/utils/writemake.tcsh > Makefile
make gradh
#
# set up the run for the first time
#
if [ $jobtype != 'restart' ]; then
cat ../setup_part1.txt > setup$x.txt
echo $x >> setup$x.txt
cat ../setup_part2.txt >> setup$x.txt
cp ../inspho .
./sph_tree_lf_gradh initial imboss < setup$x.txt
fi
#
# the following part writes a script to run the job under xgrid and submits the job
#  retaining the job id
#
rm $newdir.sh
echo '#!'/bin/tcsh > $newdir.sh;
echo cd `pwd -P` >> $newdir.sh;
echo echo hello '$HOST `date` >' hello.out >> $newdir.sh;
echo ./sph_tree_lf_gradh evolution imboss '>&' $newdir.output >> $newdir.sh
chmod a+x $newdir.sh
#
# for a dummy run do not actually submit the script
#
if [ $jobtype != 'dummy' ]; then
jobid=`xgrid $xgridauth -job submit ./$newdir.sh`
echo $jobid
jobid=${jobid/\{jobIdentifier =/};
jobid=${jobid/'; }'/};
###echo $jobid \'
#
# the following part writes scripts to clean up xgrid doo-doos.
#
echo echo checking status of $newdir... >> ../checkstatus
echo xgrid $xgridauth -job attributes -id $jobid >> ../checkstatus
echo echo getting results of $newdir... >> ../getresults
echo xgrid $xgridauth -out $PWD -se $PWD/$newdir.xgriderr -so $PWD/$newdir.xgridout -job results -id $jobid >> getresults$newdir
echo source $newdir/getresults$newdir >> ../getresults
echo echo deleting $newdir... >> ../cleanup
echo xgrid $xgridauth -job delete -id $jobid -gid $gid >> ../cleanup
#elif [ $jobtype == 'ssh' ]; then
#   echo $#machines;
#   for machine in ${machines[@]}; do
#       loadav=`ssh $machine uptime | cut -f5 -d':' | cut -f1 -d','`;
#       echo $machine load average = $loadav;
#       ssh $machine 'cd $PWD/$newdir nice +19 ./$newdir.sh < /dev/null >& crap.out &';
#   done
fi
cd ..
done;

fi
