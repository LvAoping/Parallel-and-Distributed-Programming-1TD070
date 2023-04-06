###############################################################################
# Makefile for assignment 1, Parallel and Distributed Computing 2022.
###############################################################################

CC = mpicc
CFLAGS = -std=c99 -g -O3
LIBS = -lm

BIN = stencil

all: $(BIN)

stencil: stencil.c stencil.h
	$(CC) $(CFLAGS) -o $@ $< $(LIBS)

testtest: 
	mpirun --bind-to none -np 1 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 2 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 4 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 8 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 12 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 16 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 20 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 24 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 28 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2
	mpirun --bind-to none -np 32 stencil /home/maya/public/PDP_Assignment1/input100000000.txt result.txt 2

clean:
	$(RM) $(BIN)
