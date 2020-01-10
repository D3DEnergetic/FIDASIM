title: Running FIDASIM

#Running FIDASIM

[TOC]

Running FIDASIM is as easy as running

```
lstagner@computer:~/FIDASIM$ ./fidasim
   ____ ____ ___   ___    ____ ____ __  ___
  / __//  _// _ \ / _ |  / __//  _//  |/  /
 / _/ _/ / / // // __ | _\ \ _/ / / /|_/ /
/_/  /___//____//_/ |_|/___//___//_/  /_/

Version: v2.0.0-dev

FIDASIM is released as open source code under the MIT Licence.
For more information visit http://d3denergetic.github.io/FIDASIM/
```
Actually having FIDASIM produce something takes a bit of thought

##Recommended Hardware for OpenMP Parallelized FIDASIM

The following settings will give a reasonable runtime.

* 32 threads on a shared memory node (All calculations are done on a single node)
* At least 2 GB of memory (amount of memory varies depending on user settings and inputs)

@warning
By default FIDASIM will use all the threads available. If another process is hogging a core it will cause FIDASIM to stall. To prevent this use the `num_threads` optional argument

##Recommended Hardware for MPI Parallelized FIDASIM

The following settings will give a reasonable runtime.

* 32 processes
* At least 2 GB of memory per process (amount of memory varies depending on user settings and inputs)

##Running Interactively

For OpenMP parallelized FIDASIM run with 16 threads run
```
[lstagner@dawson061]% ./fidasim /p/fida/lstagner/TEST/test_inputs.dat 16
```

For MPI parallelized FIDASIM run with 16 processes run
```
[lstagner@dawson061]% mpirun -np 16 ./fidasim /p/fida/lstagner/TEST/test_inputs.dat
```

##Submitting to a cluster job schedular
The recommended way of running FIDASIM is through a job schedular such as Slurm or PBS. 

FIDASIM provides `submit_fidasim`: a python script that schedules a FIDASIM job on a cluster. For example

```bash
lstagner@computer:~$ submit_fidasim /u/lstagner/TEST
```
will submit any incomplete FIDASIM runs in the `/u/lstagner/TEST` directory. Alternatively
```bash
lstagner@computer:~$ submit_fidasim /u/lstagner/TEST/test_inputs.dat
```
will submit just the `test` FIDASIM run.
Slurm and PBS resource managers are supported. `submit_fidasim` works for both OpenMP and MPI. Run `submit_fidasim -h` for the full documentation.

