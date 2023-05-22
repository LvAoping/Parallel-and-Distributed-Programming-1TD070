INPUT = input/32x32.txt
N = 4

CC = mpicc
CFLAGS = -std=c99 -O3
LIBS = -lm

BIN = shearsort

all: $(BIN)

shearsort: shearsort.o
	$(CC) $(CFLAGS) $^ -o $@ $(LIBS) 

clean:
	$(RM) $(BIN) *.o

%.o: %.c
	$(CC) $(CFLAGS) -c $<

run: shearsort
	mpirun -n $(N) ./shearsort -cs $(INPUT)