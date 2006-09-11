#!/bin/sh
#
# script to run a sequence of runs via xgrid
#
rm cleanup
rm getresults
for x in `cat fluxtomass_values`; do 
newdir=mbossbod_f$x
mkdir $newdir
cd $newdir
echo '--------' $newdir '--------'
#
# this part writes the input file for each run
#
pwd
cat ../setup_part1.txt > setup$x.txt
echo 40515 >> setup$x.txt
cat ../setup_part2.txt >> setup$x.txt
echo $x >> setup$x.txt
cat ../setup_part3.txt >> setup$x.txt
cp ../inspho .
~/sphNG/utils/writemake.tcsh > Makefile
make gradh
./sph_tree_lf_gradh initial imboss < setup$x.txt
#
# the following part writes a script to run the job under xgrid and submits the job
#  retaining the job id
#
rm $newdir.sh
echo '#!'/bin/sh > $newdir.sh;
echo cd $PWD >> $newdir.sh;
echo echo hello '>' hello.out >> $newdir.sh;
echo ./sph_tree_lf_gradh evolution imboss '>&' $newdir.output >> $newdir.sh
chmod a+x $newdir.sh
jobid=`xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -job submit ./$newdir.sh`
echo $jobid
jobid=${jobid/\{jobIdentifier =/};
jobid=${jobid/'; }'/};
echo $jobid
#
# the following part writes scripts to clean up xgrid doo-doos.
#
echo echo getting results of $newdir... >> ../getresults
echo xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -out $PWD -se $PWD/$newdir.xgriderr -so $PWD/$newdir.xgridout -job results -id $jobid >> getresults$newdir
echo source $newdir/getresults$newdir >> ../getresults
echo echo deleting $newdir... >> ../cleanup
echo xgrid -hostname cytosine.ex.ac.uk -auth Kerberos -job delete -id $jobid >> ../cleanup
cd ..
done;
