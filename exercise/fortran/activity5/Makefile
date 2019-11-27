FC=ftn
FFLAGS=-g

all: collective
collective: 
	$(FC) $(FFLAGS) -c collective.f90 
	$(FC) collective.o -o collective
clean:
	rm collective.o collective
