CC=cc
CFLAGS=-g

all: heat_mpi
heat_mpi: 
	$(CC) $(CFLAGS) -c heat_mpi_2d.c 
	$(CC) $(CFLAGS) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm heat_mpi_2d.o heat_mpi_2d
