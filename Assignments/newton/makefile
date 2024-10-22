# Comments explaining the structure of the makefile.

# Use variables for flexibility.
CC = gcc
CFLAGS = -lm -O3  -march=native  -lpthread -ffast-math -flto 

# Default target when no target is specified.
.PHONY: all test clean
all: newton

# Target for building the newton_iteration executable.
newton: newton.c
	$(CC) -o newton newton.c $(CFLAGS)

# Target to create a tarball of the source and makefile.
test:
	tar -czvf newton.tar.gz newton.c makefile

# Target to clean up generated files.
clean:
	rm -rvf *.o newton newton.tar.gz extracted gmon.out

RmImages:
	rm -rvf newton_convergence_* newton_attractors_*
