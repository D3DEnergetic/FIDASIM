title: Troubleshooting

This page documents the common problems people encounter while running FIDASIM.
Please try these solutions before [opening an issue](https://github.com/D3DEnergetic/FIDASIM/issues/new) on Github.

[TOC]

---

# Segfaults
If you encounter a segfault make sure you have set the stacksize limit to unlimited.

To do this for the bash shell run
```bash
ulimit -s unlimited
```
or for the tcsh shell run
```tcsh
limit stacksize unlimited
```

# Execution Hangs when using OpenMP
By default FIDASIM will use all the threads available when using OpenMP.
If another process is hogging a core it will cause FIDASIM to stall.
To prevent this use the `num_threads` optional argument as shown below

```bash
fidasim ./test_inputs.dat 8
```
