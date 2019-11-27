FC=ftn
FFLAGS=-g

all: hello
hello: 
	$(FC) $(FFLAGS) -c hello.f90 
	$(FC) hello.o -o hello
clean:
	rm hello.o hello
