FC=ftn
FFLAGS=-g

all: ringcart
ringcart: 
	$(FC) $(FFLAGS) -c ringcart.f90 
	$(FC) ringcart.o -o ringcart
clean:
	rm ringcart.o ringcart
