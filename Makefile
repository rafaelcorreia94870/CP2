CC = g++
SRC = src/
CFLAGS = -g -Wall -pg -O2 -ftree-vectorize -msse4 -mfpmath=sse -march=x86-64 -fno-omit-frame-pointer

.DEFAULT_GOAL = all

MD.exe: $(SRC)MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

perf:
	sbatch scripts/perf.sh

gprof:
	sbatch scripts/gprof.sh

a:
	sbatch scripts/NaoFUNFA.sh

stor:
	sbatch scripts/OLAATODOS.sh

clean:
	find . -type f \( ! -path "./src/*" ! -path "./scripts/*" ! -name "Makefile" ! -name "inputdata.txt" ! -name "cp_output_original.txt" ! -name "cp_output_5000.txt" ! -name "cp_average_5000.txt" ! -name "perf.sh" ! -name "run.sh" ! -name "gprof.sh" ! -name "cp_average_original.txt" ! -name "CP_Relatório.pdf" \) -exec rm -v {} \;

run:
	sbatch scripts/run.sh

run2:
	sbatch scripts/run2.sh


all: MDseq.exe MDpar.exe

MDseq.exe: $(SRC)/MDseq.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDseq.cpp -lm -o MDseq.exe

MDpar.exe: $(SRC)/MDpar.cpp
	module load gcc/11.2.0;\
	$(CC) $(CFLAGS) $(SRC)MDpar.cpp -lm -fopenmp -o MDpar.exe

runseq: MDseq.exe
	./MDseq.exe < inputdata.txt

timeseq: MDseq.exe
	./MDseq.exe < inputdata.txt

runpar: MDpar.exe
	export OMP_NUM_THREADS=2;\
	./MDpar.exe < inputdata.txt

timepar: MDpar.exe
	export OMP_NUM_THREADS=40;\
	time ./MDpar.exe < inputdata.txt