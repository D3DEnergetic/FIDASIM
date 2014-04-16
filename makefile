#System:   64 bit
#Compiler: gfortran or ifort

ifeq ($(FIDASIM_COMPILER),gfortran)
	LFLAGS = -lnetcdff -lnetcdf -lm
	CFLAGS = -g -O2 -fopenmp -Wall
endif

ifeq ($(FIDASIM_COMPILER),ifort)
	LFLAGS = -lnetcdff -lnetcdf -limf -lm
	CFLAGS = -O2 -openmp -openmp-report -parallel
endif

fidasim: fidasim.o
	$(FIDASIM_COMPILER) $(CFLAGS) fidasim.o -o fidasim -L$(NETCDF_LIB) $(LFLAGS)

fidasim.o: fidasim.f90
	$(FIDASIM_COMPILER) $(CFLAGS) -c -I$(NETCDF_INCLUDE) fidasim.f90

clean:
	-rm application.mod fidasim.o fidasim

