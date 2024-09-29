.PHONY: all test clean
all: distances


distances : distances.c
	gcc -o distances distances.c -fopenmp -O3 -lm -ffast-math -march=native -funroll-loops

test:
	tar -czvf cell_distances.tar.gz distances.c makefile
clean:
	rm -rvf *.o distances distances/ extracted/ cell_distances.tar.gz
