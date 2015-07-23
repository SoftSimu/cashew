Please refer to the manual "cashew_manual032.pdf" for detailed instructions on how to set up the program.

The program can be compiled using a command such as

f90 -o cashew_serial Mersenne.F90 Params.F90 Parallel.F90 Quaternions.F90 Functions.F90 Model.F90 IO.F90 Cashew.F90

where "f90" is the command invoking the fortran compiler. This creates a serial executable. The MPI-parallel version can be compiled with

f90 -D MPI -o cashew_parallel Mersenne.F90 Params.F90 Parallel.F90 Quaternions.F90 Functions.F90 Model.F90 IO.F90 Cashew.F90

A template Makefile is also provided for compiling using the make command.

The original source code is contained in the folder "source". The program also uses an external random number generator (RNG), which is contained in the folder "external". This is the RNG "Mersenne Twister" and it is not original source code by the authors, but obtained from http://www.math.sci.hiroshima-u.ac.jp/∼m-mat/MT/VERSIONS/FORTRAN/mt95.f90 . Here we have included it for the convenience of the librarian, editors and referees. For compiling, please move the file Mersenne.F90 in the same folder as the rest of the source code.

The folder "tutorial" contains three examples of running the program. Please see "tutorial.pdf" for instructions.

The code has a permanent web address at https://github.com/thynnine/cashew .
