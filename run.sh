#!/bin/bash
#SBATCH --time=3:00
#SBATCH --partition=cpar
#SBATCH --constraint=k20
nvcc --version
nvprof ./bin/MDcuda