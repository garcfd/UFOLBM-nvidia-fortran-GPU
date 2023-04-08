# nvidia-fortran-GPU-LBM-solver
nvidia pgf77 fortran openacc gpu solver

The laplace06 program was the beginning of getting something to run on GPU,
The main program you want is the lbm06, which is the Lattice Bolzmann solver.
The compiler commands are included in compile file, which makes 3 versons.
You will obviously need a Telsa GPU if you want to run the GPU version,
  but the other 2 versions (basic and multicore) will run fine on CPU
There are also some instructions on what you need to download and install 
to get nvidia (free of charge) installed on your machine. 

The mesh is a uniform cartesian mesh with cell size=1, so you just need 
to supply the i,j,k coordinates where the obstacle (obs) array =1.
Its read from obstacle.dat. The program contains an example of a sphere 
obstacle which can be used instead if you want to uncomment that part.

There is also a sample obstacle.dat file here with size (500 100 200)
So if you change the values of (ni nj nk) in the code it should work.
https://drive.google.com/file/d/1nfoIAPawoUsyy-OZ9c6MpxROG2A4bCsC

