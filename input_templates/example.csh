#!/bin/csh
#PBS -q normal
#PBS -l vmem=01GB
#PBS -l walltime=1:00:00


$1 = ff #filename

ifort -c -g $1.f90 >& compiler.txt
if ( $status != 0 ) then
  echo "errors compiling $1.f90"
  exit
endif
rm compiler.txt
#
ifort $1.o -lnetcdff -lnetcdf
if ( $status != 0 ) then
  echo "errors linking and loading $1.o"
  exit
endif
rm $1.o
#
mv a.out $1
./$1 > $1.txt
if ( $status != 0 ) then
  echo "errors running program"
  exit
endif
rm $1
#
echo "results written to $1.txt"
