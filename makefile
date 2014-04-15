#System:   amd64_sles11
#Compiler: ifort  ('ifort -help > man.dat' to store help into file)
COMPILER=ifort
#COMPILER=gfortran

ifeq ($(COMPILER),gfortran)
	LFLAGS = -lnetcdff -lnetcdf -lm
	CFLAGS = -g -O2 -fopenmp -Wall -fbacktrace
endif

ifeq ($(COMPILER),ifort)
	LFLAGS = -lnetcdf -lnetcdf -limf -lm
	CFLAGS = -O2 -openmp -openmp-report -parallel
endif

fidasim: fidasim.o
	$(COMPILER) $(CFLAGS) fidasim.o -o fidasim -L$(NETCDF_LIB) $(LFLAGS)

fidasim.o: fidasim.f90
	$(COMPILER) $(CFLAGS) -c -I$(NETCDF_INCLUDE) fidasim.f90

clean:
	-rm application.mod fidasim.o fidasim

