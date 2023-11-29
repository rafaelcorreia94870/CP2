#!/bin/bash
#SBATCH --time=10:00
#SBATCH --cpus-per-task=40
#SBATCH --partition=cpar
perf stat ./MDseq.exe < inputdata.txt
perf stat ./MDpar.exe < inputdata.txt 
