FC=gfortran 
FFLAGS=-O3 -Wall -Wextra
#FFLAGS = -g  -fbounds-check -Waliasing -Wall -fbacktrace -ffpe-trap=invalid,overflow,zero  -finit-real=nan -finit-integer=nan  -ffree-line-length-none  -O0 -Warray-bounds  -ffast-math -funroll-loops -ftree-vectorizer-verbose=2
#FFLAGS = 
LOADFLAG= -static

#EXEC=../exec/scm_pgmet.exe
EXEC=./scm_pgmet.exe

OBJ=  Main.o  Class_Module_Constants.o Class_Module_Fields.o Class_Module_Dynamics.o \
              Class_Module_TimeManager.o

Main:	$(OBJ)
	$(FC) $(FFLAGS) -o $(EXEC) $(LOADFLAG) $(OBJ)
  
Class_Module_Constants.o: Class_Module_Constants.f90
	$(FC) $(FFLAGS) -c Class_Module_Constants.f90

Class_Module_TimeManager.o: Class_Module_TimeManager.f90 Class_Module_Constants.o
	$(FC) $(FFLAGS) -c Class_Module_TimeManager.f90

Class_Module_Fields.o: Class_Module_Fields.f90 Class_Module_Constants.o
	$(FC) $(FFLAGS) -c Class_Module_Fields.f90

Class_Module_Dynamics.o: Class_Module_Dynamics.f90 Class_Module_Fields.o Class_Module_Constants.o
	$(FC) $(FFLAGS) -c Class_Module_Dynamics.f90

Main.o: Main.f90  Class_Module_Constants.o Class_Module_TimeManager.o Class_Module_Dynamics.o \
	Class_Module_Fields.o
	$(FC) $(FFLAGS) -c Main.f90

.SUFFIXES:
.SUFFIXES: .f90 .o

.f90.o:
	$(FC) $(FFLAGS) -c $@ $<

clean:
	-rm -f $(OBJ)
	-rm -f $(EXEC)
	-rm -f *.mod
	-rm -f *.o
