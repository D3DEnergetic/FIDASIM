#System:   amd64_sles11
#Compiler: ifort  ('ifort -help > man.dat' to store help into file)
COMPILER=ifort -O2 -openmp -openmp-report -parallel
NETCDF_INCLUDE=/usr/pppl/intel/11-pkgs/netcdf-4.1.3/include
NETCDF_LIB=/usr/pppl/intel/11-pkgs/netcdf-4.1.3/lib
#COMPILER=gfortran44 -fopenmp -Wall 
#COMPILER=pgf90 -O2 -Mconcur -Minform=inform 

fidasim: fidasim.o
	$(COMPILER) fidasim.o -o fidasim -L$(NETCDF_LIB) -lnetcdff -lnetcdf -limf -lm

fidasim.o: fidasim.f90
	$(COMPILER) -c -I$(NETCDF_INCLUDE) fidasim.f90

clean:
	-rm application.mod fidasim.o fidasim

