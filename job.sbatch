#!/bin/bash 

#SBATCH --nodes=3                        # requests 3 compute servers
// #SBATCH --ntasks-per-node=1              # runs 2 tasks on each server
#SBATCH --cpus-per-task=1                # uses 1 compute core per task
#SBATCH --time=1:00:00
#SBATCH --mem=2GB
#SBATCH --job-name=int_ring.cpp
#SBATCH --output=int_ring.out

module purge
module load  amber/openmpi/intel/20.06 

mpiexec ./int_ring