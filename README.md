# nvidia-fortran-GPU-LBM-solver
nvidia pgf77 fortran openacc gpu solver

This code was developed from the 2D python code of Jonas Latt (unige coursera).
Mainly converting it to (D3Q27) in Fortran and making it run on nvidia GPU.
The compiler commands are included in "compile" file, giving 3 executables.

You will obviously need a Telsa GPU if you want to run the GPU version,
but the other 2 versions (basic and multicore) will run fine on just CPU.
There are also some instructions on what you need to download and install 
to get nvidia (free of charge) installed on your machine. 

The mesh is a uniform cartesian mesh with cell size=1, so you just need 
to supply the i,j,k coordinates where the obstacle (obs) array=1. There is 
a sample obstacle.dat file for reference. The program contains an example of 
a sphere obstacle which can be used instead if you want to uncomment that part.

Recent work has focussed on the D3Q27 cell type, and developing
a reflection matrix which can be used to model free-slip walls.
(https://www.youtube.com/watch?v=Y9qtjCZnT-U)
