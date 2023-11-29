#!/bin/bash
#SBATCH --partition=cpar
#SBATCH --cpus-per-task=40
#SBATCH --time=10:00
for cores in 1 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40; do
    time srun --partition=cpar --cpus-per-task=$cores ./MDpar.exe < inputdata.txt 
done