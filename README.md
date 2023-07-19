# nvidia-fortran-GPU-LBM-solver
nvidia pgf77 fortran openacc gpu solver

The laplace-06 program was just the start of getting something to run on GPU,
The main program you want to use is the lbm-06, which is the Lattice Bolzmann solver.
The compiler commands are included in "compile" file, which makes 3 executable files.
You will obviously need a Telsa GPU if you want to run the GPU version,
  but the other 2 versions (basic and multicore) will run fine on just CPU.
There are also some instructions on what you need to download and install 
to get nvidia (free of charge) installed on your machine. 

The mesh is a uniform cartesian mesh with cell size=1, so you just need 
to supply the i,j,k coordinates where the obstacle (obs) array=1.
Its read from obstacle.dat. The program contains an example of a sphere 
obstacle which can be used instead if you want to uncomment that part.

There is also a sample obstacle.dat file here with size (500 100 200).
So if you change the values of (ni nj nk) in the code it should work.
https://drive.google.com/file/d/1nfoIAPawoUsyy-OZ9c6MpxROG2A4bCsC

[There is also porsche.stl file used to create the obstacle file.
https://drive.google.com/file/d/1_Z8FRGjBYJTyB6_KUXJ2Iy8ZvH4JUBLG

Recent work has focussed on the D3Q27 cell type, and developing
a reflection matrix which can be used to model free-slip walls.
(https://www.youtube.com/watch?v=Y9qtjCZnT-U)
