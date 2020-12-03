#
# Makefile for ct-cbn
#

VERSION = 0.1.00

#CC = gcc-4.2 -Wall -O3 -march=core2 -fopenmp  #-pg
CC = gcc -Wall -O3 -fopenmp -pg -mcmodel=medium -lbsd
#CC = gcc-9 -Wall -O3 -fopenmp# -pg


all: clonex

clonex.o: clonex.c
	$(CC)  -I/usr/local/include -c clonex.c -lbsd

clonex: clonex.o
	$(CC) clonex.o -L/usr/local/lib -lgsl -lgslcblas -lm -lbsd -o $@

clean:
	rm -f a.out core *.o clonex *~ \#*



