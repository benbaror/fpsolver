#F90=ifort
#F90=g95
F90=ifort
#FLAGS = -O3
FLAGS=-i-dynamic -O3 -g -pg #-mpe=log

#
# Set preprocessor variable to include loss-cone term 
#
LCflag = 
LCflag= -D__LCTERM__
PPA = # Preprocessor arguments. Usege: make "PAA = -D (macro)"
all: FPSolverMS #move_mod


FPSolverMS: ./temp/FPSolverMS.o ./temp/FPglobal.o ./temp/FPvars.o ./temp/FPinit.o ./temp/FPgrid.o ./temp/FPwork.o ./temp/FPcurrent.o ./temp/FPpdecol.o ./temp/FPuser.o
	 $(F90) $(FLAGS) -o FPSolverMS  ./temp/FPSolverMS.o ./temp/FPglobal.o ./temp/FPvars.o ./temp/FPinit.o ./temp/FPgrid.o ./temp/FPwork.o ./temp/FPcurrent.o ./temp/FPpdecol.o ./temp/FPuser.o 

#move_mod:
#	mv *.mod ./temp
# Here are the compile steps

#
# Version w/o LC
#
./temp/FPSolverMS.o:FPSolverMS.F90 ./temp/FPglobal.o ./temp/FPvars.o ./temp/FPinit.o ./temp/FPwork.o ./temp/FPgrid.o ./temp/FPcurrent.o ./temp/FPpdecol.o
	 $(F90) $(FLAGS) $(LCflag) $(PPA) -c  FPSolverMS.F90 -o ./temp/FPSolverMS.o

./temp/FPglobal.o:./FPglobal.f90
	 $(F90) $(FLAGS) -c  ./FPglobal.f90 -o ./temp/FPglobal.o

./temp/FPvars.o:./FPvars.f90
	 $(F90) $(FLAGS) -c  ./FPvars.f90 -o ./temp/FPvars.o

./temp/FPinit.o:./FPinit.f90 ./temp/FPglobal.o ./temp/FPvars.o ./temp/FPgrid.o
	 $(F90) $(FLAGS) -c  ./FPinit.f90 -o ./temp/FPinit.o

./temp/FPgrid.o:./FPgrid.f90
	 $(F90) $(FLAGS) -c  ./FPgrid.f90 -o ./temp/FPgrid.o

./temp/FPwork.o:./FPwork.f90
	 $(F90) $(FLAGS) -c  ./FPwork.f90 -o ./temp/FPwork.o

./temp/FPcurrent.o:./FPcurrent.f90 ./temp/FPglobal.o ./temp/FPvars.o
	 $(F90) $(FLAGS) -c  ./FPcurrent.f90 -o ./temp/FPcurrent.o

./temp/FPpdecol.o:./FPpdecol.f90 ./temp/FPglobal.o ./temp/FPwork.o ./temp/FPuser.o ./temp/FPgrid.o
	 $(F90) $(FLAGS) -c  ./FPpdecol.f90 -o ./temp/FPpdecol.o

./temp/FPuser.o:./FPuser.F90 ./temp/FPglobal.o ./temp/FPvars.o ./temp/FPgrid.o ./temp/FPcurrent.o
	 $(F90) $(FLAGS)  $(LCflag) $(PPA) -c  ./FPuser.F90 -o ./temp/FPuser.o

# This entry allows you to type " make clean " to get rid of
# all object and module files
clean:
	 rm -f -r f_{files,modd}* *.o *.mod *.M *.d V*.inc *.vo  V*.f *.dbg album F.err temp/*.*
