#!/bin/bash



#SBATCH --job-name=int_ring
#SBATCH --nodes=6 
#SBATCH --cpus-per-task=1
#SBATCH --mem=2GB 
#SBATCH --time=01:00:00 
#SBATCH --verbose

module purge 
module load  amber/openmpi/intel/20.06 


mpiexec ./int_ring
