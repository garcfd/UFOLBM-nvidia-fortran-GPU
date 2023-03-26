# nvidia-fortran-GPU-LBM-solver
nvidia pgf77 fortran openacc gpu solver

The laplace06 program was the beginning of getting something to run on GPU,
The main program you want is the lbm06, which is the Lattice Bolzmann solver.
The compiler commands are included in compile file, which makes 3 versons.
You will obviously need a Telsa GPU if you want to run the GPU version,
  but the other 2 versions (basic and multicore) will run fine on CPU
There are also some instructions on what you need to download and install 
to get nvidia (free of charge) installed on your machine. 
Good luck and ask questions if you need to (garcfd@gmail.com)
