CPUCC       = mpic++ #/opt/gnu-gcc/bin/g++
CPPFLAGS    = --std=c++11 -m64 -Wno-reorder #     -Wunused-variable #-Wno-conversion-null -Wdelete-non-virtual-dtor
DEFS        = $(APPLICATION)  $(VERSION) -DALIGN=64

INCLUDES  =
LIBDIR   =
#LIBDIR      = -L/opt/gnu-gcc/lib64
LIBS     =

LINK     =  $(LIBDIR) $(LIBS)
CPU_COMPILE = $(CPUCC) $(DEFS) $(INCLUDES) $(CPPFLAGS)


all:
	$(CPU_COMPILE) createMatrix.cpp $(LINK) -o createMatrix
	$(CPU_COMPILE) loadMatrix.cpp $(LINK) -o loadMatrix
	$(CPU_COMPILE) transpose.cpp $(LINK) -o transpose

transposeMatrix:
	rm -f transpose mat-out.bin
	$(CPU_COMPILE) transpose.cpp $(LINK) -o transpose
	mpirun -np 2 ./transpose

clean:
	rm -f createMatrix loadMatrix transpose
