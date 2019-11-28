#FC = gfortran
FC = ifort 
#FCFLAGS = -g -Wall -fdefault-real-8
FCFLAGS = -qopenmp
#INCLUDES = -I/home/wktang/hdf5-1.8.15/include
#LFLAGS = -L/home/wktang/hdf5-1.8.15/lib
#LIBS= -lhdf5_fortran
SO =  tm_slab.o var.o initial.o boundary.o \
      plot1.o exitinfo.o smooth.o rhs.o rk4.o \

main: $(SO)
	$(FC) $(FCFLAGS) $(LFLAGS) $(LIBS) $(SO) -o main

%.o: %.f90 var.o
	$(FC) $(FCFLAGS) $(INCLUDES) -c $<

var.o:var.f90
	$(FC) $(FCFLAGS) $(INCLUDES) -c $<
	
clean:
	rm -f errfile outfile *.TXT *.o  
