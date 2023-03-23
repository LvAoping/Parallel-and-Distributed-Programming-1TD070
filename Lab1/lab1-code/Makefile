############################
# Makefile for lab 1
# Author: Malin Kallen 2019
############################

CC = mpicc
CCFLAGS = -g -std=c99 -O3
LIBS = -lm

BINS = hello deadlock exchange onetoall pingpong pi buffon-laplace

all: $(BINS)

buffon-laplace: buffon-laplace.c buffon-laplace.h
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

%: %.c
	$(CC) $(CCFLAGS) -o $@ $< $(LIBS)

clean:
	$(RM) $(BINS)

