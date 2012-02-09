#!/bin/csh
#PBS -l ncpus=16
#ncpus must be a multiple of 16
#PBS -l walltime=5:00:00
#PBS -q batch
#PBS -N Conv-160-4
#PBS -M leftynm@umich.edu
#PBS -m ea
# -m email when: a=abort b=begin e=end s=suspend n=none
#PBS -j oe
#PBS -o $HOME/jobs
#PBS -e $HOME/jobs

set echo

ja

#cd to directory qsub was called from
cd $PBS_O_WORKDIR

#run my executable
##ulimit -s unlimited
setenv OMP_NUM_THREADS $PBS_NCPUS
#setenv KMP_STACKSIZE 4G
omplace -nt $OMP_NUM_THREADS ./expl.out

ja -chlst