# FIDASIM
FIDASIM is a code that models the signal that is produced by charge-exchange reactions between fast-ions and injected neutral beams in tokamak plasmas. 

# How to Install 
Warning! The following procedure creates a directory named FIDASIM. This will overwrite any directory also named FIDASIM. Take necessary precautions. 

## DIII-D
On venus

    git clone https://github.com/D3DEnergetic/FIDASIM.git FIDASIM
    cd FIDASIM

To compile the stable version run 

    git checkout master
    make -f makefile_d3d

else to compile the development version run

    git checkout development
    make -f makefile_d3d

## NSTX-U
On portal

    module load git
    module load intel
    git clone https://github.com/D3DEnergetic/FIDASIM.git FIDASIM
    cd FIDASIM

To compile the stable version run 

    git checkout master
    make -f makefile_nstx

else to compile the development version run

    git checkout development
    make -f makefile_nstx
	
## References

Heidbrink, W. W., et al. "A code that simulates fast-ion D-alpha and neutral particle measurements." Comm. Comp. Phys. 10 (2011) 716.

Geiger, Benedikt. "Fast-ion transport studies using FIDA spectroscopy at the ASDEX Upgrade tokamak." Diss. lmu, 2013. APA	

