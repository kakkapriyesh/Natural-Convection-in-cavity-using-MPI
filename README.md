# kakkapriyesh-MPI_in_Python
These codes will be useful for engineers (CFD/FEM) to learn MPI.The code is universal for any aspect ratio and other parameters. MPI does not support Python and the underlying library for MPI is in C++.
Following are the codes:
 a) Stream_vorticity_formulations: Heated cavity is being solved using stream vorticity formulation and parallel processors.
 b) Tag, Ring are the methods to pass data among the processors.
 c) Matrix code can be used to solve large Matrix with parallel processors. With 4 procs I have got a speed up of nearly 3x.
 d) Coutte flow can be used to understand parallization of 2d jocobi and how boundary condition is delt with.
