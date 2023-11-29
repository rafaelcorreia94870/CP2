#!/bin/bash
#SBATCH --time=10:00
#SBATCH --cpus-per-task=40
#SBATCH --partition=cpar
./MDpar.exe < inputdata.txt
gprof MDpar.exe gmon.out > analysis1.txt
cat analysis1.txt