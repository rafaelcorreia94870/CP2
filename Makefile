CC = g++
SRC = src/
CFLAGS = -Wall -pg -O2 -ftree-vectorize -mavx -mfpmath=sse -march=x86-64

.DEFAULT_GOAL = MD.exe

MD.exe: $(SRC)MD.cpp
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe

test:
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe
	srun --partition=cpar perf stat -M cpi ./MD.exe < inputdata.txt

testrepeat:
	$(CC) $(CFLAGS) $(SRC)MD3.cpp -lm -o MD3.exe
	srun --partition=cpar perf stat -r 10 -e instructions,cycles ./MD3.exe < inputdata.txt

gprof:
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe
	./MD.exe < inputdata.txt
	gprof MD.exe gmon.out > analysis1.txt
	cat analysis1.txt

clean:
	find . -type f \( ! -path "./src/*" ! -name "Makefile" ! -name "inputdata.txt" ! -name "cp_output_original.txt" ! -name "cp_average_original.txt" ! -name "CP_Relatório.pdf" \) -exec rm -v {} \;

run:
	srun  --partition=cpar ./MD.exe < inputdata.txt
