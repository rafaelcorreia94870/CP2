CC = g++
SRC = src/
CFLAGS = -Wall -pg -O2 -ftree-vectorize -mavx -mfpmath=sse -march=x86-64

.DEFAULT_GOAL = MD.exe

test:
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe
	srun --partition=cpar perf stat -M cpi ./MD.exe < inputdata.txt

test2:
	$(CC) $(CFLAGS) $(SRC)MD2.cpp -lm -o MD2.exe
	srun --partition=cpar perf stat -M cpi ./MD2.exe < inputdata.txt

test3:
	$(CC) $(CFLAGS) $(SRC)MD3.cpp -lm -o MD3.exe
	srun --partition=cpar perf stat -M cpi ./MD3.exe < inputdata.txt

testrepeat:
	$(CC) $(CFLAGS) $(SRC)MD3.cpp -lm -o MD3.exe
	srun --partition=cpar perf stat -r 3 -M cpi ./MD3.exe < inputdata.txt

gprof1:
	$(CC) $(CFLAGS) $(SRC)MD.cpp -lm -o MD.exe
	./MD.exe < inputdata.txt
	gprof MD.exe gmon.out > analysis1.txt
	cat analysis1.txt

gprof2:
	$(CC) $(CFLAGS) $(SRC)MD2.cpp -lm -o MD2.exe
	./MD2.exe < inputdata.txt
	gprof MD2.exe gmon.out > analysis2.txt
	cat analysis2.txt

gprof3:
	$(CC) $(CFLAGS) $(SRC)MD3.cpp -lm -o MD3.exe
	./MD3.exe < inputdata.txt
	gprof MD3.exe gmon.out > analysis3.txt
	cat analysis3.txt

clean:
	find . -type f \( ! -path "./src/*" ! -name "Makefile" ! -name "inputdata.txt" ! -name "cp_output_original.txt" ! -name "cp_average_original.txt" \) -exec rm -v {} \;

run:
	./MD.exe < inputdata.txt
