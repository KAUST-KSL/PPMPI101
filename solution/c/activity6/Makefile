CC=cc
CFLAGS=-g

all: heat_serial heat_mpi

heat_serial: 
	$(CC) $(CFLAGS) -c heat_serial.c 
	$(CC) $(CFLAGS) heat_serial.o -o heat_serial

heat_mpi: 
	$(CC) $(CFLAGS) -c heat_mpi_2d.c 
	$(CC) $(CFLAGS) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm *.o heat_mpi_2d heat_serial
