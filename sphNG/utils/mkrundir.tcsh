#!/bin/tcsh
#
# @(#) Copies sph files to new directory and creates script
# @(#)  -- Daniel Price 11/9/03
#
if $# != 3 then
   echo Usage: $0 dirname ifile initialdumpname
   echo e.g. $0 crap icrap crapp01
else
   set newdir=$1
   echo making directory $newdir
   mkdir $newdir
   echo moving files...
   cp $2 $newdir
   cp $3 $newdir   
   cp sph_tree $newdir
   cp inspho $newdir
   cd $newdir
   ##../write_sgescript $2 > $2.sge 
   cp $2 $2\.s
## write a Makefile
   ../scripts/writemake.tcsh > Makefile
## append a rerun option to Makefile
   echo rerun: >> Makefile
   echo '	cp $2\.s $2' >> Makefile 
## link for supersphplot
   ##ln -s ../supersphplot supersphplot
   echo "that's all folks"
endif
