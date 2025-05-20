Parallel Gravitational Field Computation

This project implements a parallel algorithm to compute the gravitational field using the MPI (Message Passing Interface) library.
The code solves a system of linear equations to calculate gravitational potential and exports the results to data files.
Features

Parallel Data Loading: 
Process 0 (P0) loads input data (B, L, H, g values) from a file.
MPI Broadcast: Distributes B, L, H, and g values to all processes.
Coordinate Computation: Calculates X, S, and normal vectors (n) for each point.
Matrix Distribution: Assigns parts of the system matrix A to each process.
BSGS Solver: Implements a parallel Bi-Conjugate Gradient Stabilized (BCGS) solver.
Computes matrix-vector products in parallel.
Process 0 prints residual norms.

Data Gathering: Uses MPI_Gather or MPI_Allgather to collect alpha vectors from all processes.
Vector U Computation: Calculates the solution vector U in parallel.
Data Export: Process 0 writes alpha vectors, vector U, and final data (B, L, U) to output files.

Dependencies

MPI (e.g., MPICH or OpenMPI)
C++ Compiler (e.g., g++)
Standard math library (<cmath>)

Usage

Compile the code with MPI support:mpicxx -o TiazovePoleParalelne TiazovePoleParalelne.cpp

Run the parallel program:mpirun -np <number_of_processes> ./TiazovePoleParalelne

Ensure the input file BL-8102.dat and its location.

Output files

outputAlpha.dat: Contains alpha values from the BSGS solver.
outputVectorU.dat: Contains computed vector U values.
outputData.dat: Contains B, L, and U values for each point.

Notes

The program assumes a fixed number of points (N = 8102) and constants (R = 6378.0, D = 300, GM = 398600.5). Adjust the input file path and tolerance (TOL) as needed.
