#!/bin/bash
#SBATCH --time=10:00
#SBATCH --cpus-per-task=40
#SBATCH --partition=cpar
./MD.exe < inputdata.txt
gprof MD.exe gmon.out > analysis1.txt
cat analysis1.txt