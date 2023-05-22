CC = mpicc
CCFLAGS = -g -O3
LIBS = -lm

BINS = matmul 

all: $(BINS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

testtest:
	mpirun -np 4 ./matmul /proj/uppmax2023-2-13/nobackup/matmul_indata/input1800.txt 111.txt

clean:
	$(RM) $(BINS)