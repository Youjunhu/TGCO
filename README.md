# TGCO
Monte-Carlo simulation of neutral beam injection and fusion alpha particle heating.

Compile:
mpif90 src_in_one_file.f90 -llapack -lblas

Run:
mpirun -n 4 ./a.out     

Input: input.nmlt, gfile/g101473.04500, gfile/ne.dat, gfile/te.dat, gfile/ti.dat

Output: various .txt files
