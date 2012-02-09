#!/bin/sh

###############################
# PBS settings Here

# Which shell to use
#PBS -S /bin/sh

# Email address
#PBS -M leftynm@umich.edu

# When to send an email: a=abort b=begin e=end s=suspend n=none
#PBS -m ea

# Name of the run
#PBS -N Mixed-1-080

# Redirect the STDOUT and STDERR files to the ~/jobs directory
#PBS -o /nobackup/leftynm/jobs
#PBS -e /nobackup/leftynm/jobs
#PBS -j oe

# Which queue - change from route to cac for debugging
#PBS -q route

#PBS -A thornton

# Number of nodes (not processors) and wallclock estimate hh:mm:ss
#PBS -l nodes=1:ppn=4,feature='opt2218|opt2220',walltime=3:00:00,qos=thornton

# End PBS Settings
###############################
ulimit -c 0

myDir=$PBS_O_WORKDIR

filename=expl.out

export OMP_NUM_THREADS=4

cd $myDir/

./$filename
