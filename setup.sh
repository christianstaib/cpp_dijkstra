module load compiler/gnu
module load mpi/openmpi
export OMP_NUM_THREADS=80
export MKL_NUM_THREADS=80

# OMP_PROC_BIND=true OMP_NUM_THREADS=40 mpirun --map-by ppr:1:package main
