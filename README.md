# Threaded Orthotree
Threaded Fully-Threaded-Tree Implementation

This project aims to better distribute an adaptive mesh across processors, by using a Hilbert space filling curve. This is stored as a thread through the leaves of the mesh allowing each process to perform operations on its own section of the tree without knowledge of the rest of the structure, without complex management of the data.

# Building and Running

The project is comprised of the lib and tests folders. There are a number of tests but the one that currently best shows the functionality of the library is 'parpois' (parallel poisson). 

```
$ mkdir build
$ cmake ..
$ make parpois
$ mpirun -n 4 ./test/parpois/parpois poisson.par
```
# Parameters

In the above sequence of instructions, poisson.par is used to pass instructions to the test. For example the following is my test configuration:

```
minDepth=5
maxDepth=7
tftt.two2one=3
tftt.ghosts=1
dirichlet=0.2
plotEvery=2000
iterations=10000
```

These parameters will remain undocumented as they are highly dependant on the test and likely to change in the near future. 

Parameters can also be specified on the command line using the syntax: "+minDepth 5" for example. It is possible to specify both a file and individual parameters, in which case they will be read in the order specified. 
