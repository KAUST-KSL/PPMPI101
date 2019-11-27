#!/bin/bash -l

#SBATCH --job-name=mpi_training
#SBATCH --time=00:10:00
#SBATCH --reservation=training
#SBATCH --nodes=1

srun -n 2 ./hello

