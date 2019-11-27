FC=ftn
FFLAGS=-g

all: pingpong
pingpong: 
	$(FC) $(FFLAGS) -c pingpong.f90 
	$(FC) pingpong.o -o pingpong
clean:
	rm pingpong.o pingpong
