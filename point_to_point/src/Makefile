FC = mpiifort
FFLAGS = -O3 -fp-model precise
 
 NC_LIB = $(shell nf-config --flibs)
 NC_INC = -I$(shell nf-config --includedir)
 LDFLAGS = $(NC_INC)
 LIBS = $(NC_LIB)
 
.SUFFIXES:
.SUFFIXES: .F90 .o
 
.F90.o:
	$(FC) -c $(FFLAGS) $(LDFLAGS) $<
 
include depend.mk

clean :
	rm -f *.o *.mod *.exe 

depend depend.mk:
	makedepf90 -o wave.exe *.F90 > depend.mk
