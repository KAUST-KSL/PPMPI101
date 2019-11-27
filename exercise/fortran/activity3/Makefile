FC=ftn
FFLAGS=-g

all: ring
ring: 
	$(FC) $(FFLAGS) -c ring.f90 
	$(FC) ring.o -o ring
clean:
	rm ring.o ring
