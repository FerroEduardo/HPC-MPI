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
	$(CPU_COMPILE) 02-PingPong.cpp $(LINK) -o PingPong
	$(CPU_COMPILE) 02-PingPongI.cpp $(LINK) -o PingPongI

clean:
	rm PingPong PingPongI
