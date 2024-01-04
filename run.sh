#!/bin/bash
#SBATCH --time=3:00
#SBATCH --partition=cpar
#SBATCH --constraint=k20
time nvprof ./bin/MDcuda < inputdata.txt
