CC = g++
SRC = src/
CFLAGS = -Wall -pg -O2 -ftree-vectorize -mavx -mfpmath=sse -march=x86-64 -lm -fno-omit-frame-pointer -fopenmp

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

perf:
	sbatch perf.sh

gprof:
	sbatch gprof.sh

clean:
	find . -type f \( ! -path "./src/*" ! -name "Makefile" ! -name "inputdata.txt" ! -name "cp_output_original.txt" ! -name "cp_output_5000.txt" ! -name "cp_average_5000.txt" ! -name "perf.sh" ! -name "run.sh" ! -name "gprof.sh" ! -name "cp_average_original.txt" ! -name "CP_Relatório.pdf" \) -exec rm -v {} \;

run:
	sbatch run.sh
