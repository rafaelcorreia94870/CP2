#!/bin/bash
#SBATCH --ntasks=40
#SBATCH --time=00:10:00
#SBATCH --partition=cpar
#SBATCH --exclusive
#SBATCH --cpus-per-task=40

threads=(1 2 4 20 40)

for nthreads in "${threads[@]}"
do
    OMP_NUM_THREADS=${nthreads} srun ./MDpar.exe <inputdata.txt >lixo
    echo "Number of threads: ${OMP_NUM_THREADS}"
done
