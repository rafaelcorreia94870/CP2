#!/bin/bash
#SBATCH --time=10:00
#SBATCH --cpus-per-task=25
#SBATCH --partition=cpar
perf stat ./MDpar.exe < inputdata.txt 