#!/bin/bash -l


#SBATCH --reservation=training
#SBATCH --time=00:20:00
#SBATCH --job-name=ddt_session

module load arm-forge
ddt --connect
