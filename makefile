#FC = gfortran
FC = ifort 
#FCFLAGS = -g -Wall -fdefault-real-8
FCFLAGS = -qopenmp
#INCLUDES = -I/home/wktang/hdf5-1.8.15/include
#LFLAGS = -L/home/wktang/hdf5-1.8.15/lib
#LIBS= -lhdf5_fortran
SO =   tm_slab.o var.o initial.o  \
#      start_main.o start_initial.o nonlinear_enr.o \
      nonlinear_main.o  nonlinear_rhs.o \
      conjugates.o derivatives.o  nonlinear_pus.o replace.o  \
      solver.o contour_data.o random.o boundary.o

#main: main.o param.o dmotifs.o ssa.o
#    $(FC) $^ $(LFLAGS) $(LIBS) -o main 

#param.o: param.f90
#    $(FC) $(FCFLAGS) -c $<

#dmotifs.o: dmotifs.f90 param.o
#    $(FC) $(FCFLAGS) -c $<

#ssa.o: ssa.f90 dmotifs.o
#    $(FC) $(FCFLAGS) -c $<

#main.o: main.f90 param.o dmotifs.o ssa.o
#    $(FC) $(FCFLAGS) -c $(INCLUDES) $<

main: $(SO)
	$(FC) $(FCFLAGS) $(LFLAGS) $(LIBS) $(SO) -o main

%.o: %.f90 var.o
	$(FC) $(FCFLAGS) $(INCLUDES) -c $<

var.o:var.f90
	$(FC) $(FCFLAGS) $(INCLUDES) -c $<
	
clean:
	rm -f errfile outfile *.TXT *.o  
