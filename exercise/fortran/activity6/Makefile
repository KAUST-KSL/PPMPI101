FC=ftn
ifeq ($(PE_ENV),CRAY)
   FFLAGS=-g -eZ
else ifeq ($(PE_ENV),GNU)
   FFLAGS=-g -cpp
endif

all: heat_mpi
heat_mpi: 
	$(FC) $(FFLAGS) -c heat_mpi_2d.f90 
	$(FC) $(FFLAGS) heat_mpi_2d.o -o heat_mpi_2d
clean:
	rm heat_mpi_2d.o heat_mpi_2d *.mod *.i
