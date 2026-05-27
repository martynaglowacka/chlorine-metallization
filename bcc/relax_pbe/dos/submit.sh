#!/bin/bash -l 
#SBATCH -J pieski
#SBATCH -N 2
#SBATCH --ntasks-per-node=28
#SBATCH -A g87-1106
#SBATCH -p topola
#SBATCH -t 00:20:00

module load apps/fhi-aims/221103
srun aims.221103.hdf5.scalapack.mpi.x  control.in | tee aims.out
