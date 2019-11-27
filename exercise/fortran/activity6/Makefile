FC=ftn
FFLAGS=-g -eZ

all: heat_mpi
heat_mpi: 
	$(FC) $(FFLAGS) -c heat_mpi_2d.f90 
	$(FC) $(FFLAGS) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm heat_mpi_2d.o heat_mpi_2d *.mod *.i
