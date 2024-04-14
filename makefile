# YouJunHu's makefile
program_name=TGCO
#ifeq  (${FC}, /pkg/intel/composer_xe_2015.0.090/bin/intel64/ifort)
ifeq  (${which_machine}, parallel21)
 #   COMPILER=	ifort
   COMPILER=	mpiifort
#    OPTION= -O -fopenmp
    OPTION= -O3
lib2 = -llapack -lblas
else ifeq (${which_machine}, parallel01)
   COMPILER=	mpiifort
    OPTION= -O3
lib2 = -llapack -lblas
else ifeq (${which_machine},sugon_hefei)
   COMPILER=	mpiifort
   #OPTION= -O2 -traceback
OPTION= -O3 
   lapack = /usr/lib64/liblapack.so.3
else
#   COMPILER=	gfortran
   COMPILER=	mpif90
#    OPTION= -O -fopenmp 
   #OPTION= -fbounds-check -fopenmp
#   OPTION= -fbounds-check -g -fopenmp
#  OPTION= -Wimplicit-interface -fPIC -fmax-errors=1 -g -fcheck=all -fbacktrace -fimplicit-none -fbounds-check
  OPTION= -g -fcheck=all -fbacktrace -fbounds-check -fimplicit-none   -fno-realloc-lhs -std=f2008 -pedantic -fopenmp -Wconversion
#-Wextra 
# OPTION= -fcheck=all  -Wall -Wno-unused-variable -Wno-unused-dummy-argument -Wno-tabs -fbounds-check -fimplicit-none
# OPTION= -fcheck=all  -Wall -Wno-unused-variable -Wunused-dummy-argument -Wno-tabs -fbounds-check -fimplicit-none
lib2 = -llapack -lblas
endif

COMPILE=	$(COMPILER) $(OPTION) -c 

f90sources=src_in_one_file.f90
#solovev.f90
#f77sources= 

f90objs= $(f90sources:.f90=.o)
#f77objs= $(f77sources:.for=.o)

$(program_name): $(f90objs) $(f77objs)
	$(COMPILER)  $(OPTION)  $(f90objs) $(f77objs)  ${lapack} -o $(program_name) ${lib2}
%.o: %.f90 
	$(COMPILE) $< -o $@

.PHONY : clean run 
clean :
	 rm -f program_name  $(f90objs) $(f77objs) *.mod  $(program_name)

run: $(program_name)
	mpirun -np 4 ./$(program_name)
