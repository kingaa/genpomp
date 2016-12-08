INCLUDEDIRS = $(shell R CMD config --cppflags) $(shell pkg-config --cflags libRmath libR)
INCLUDEDIRS += -I. #-I../../src/models/$(usermodelDirectory)                                                                                                                                                                                                                  
CCFLAGS = -fopenmp -g -DDEBUG -mkl
CCFLAGS += $(CMD_LINE_FLAGS)
CFLAGS = -Wall -std=gnu99 -fpic -fopenmp -g -O2 
CFLAGS += $(CMD_LINE_FLAGS)
LFLAGS = $(shell R CMD config --ldflags) -lgfortran -lquadmath
LFLAGS += $(shell R CMD config LAPACK_LIBS)
LFLAGS += $(shell R CMD config BLAS_LIBS)
LFLAGS += $(shell pkg-config --libs libRmath)
LFLAGS += $(shell pkg-config --libs libR)
#GXX = g++
#GCC =  gcc
CXX = icpc -mkl #-fp-model precise #-fp-trap="invalid,divzero,underflow,overflow"                                                                                                                                                                                             
CC = icc -mkl #-fp-model precise #-fp-trap="invalid,divzero,underflow,overflow"                                                                                                                                                                                               

default: $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.cc,%.o,$(wildcard *.cc))

%.exe: %.cc 
	make $(patsubst %.c,%.o,$(wildcard *.c))
	$(CXX) -o $*.exe $*.cc $(patsubst %.c,%.o,$(wildcard *.c)) $(INCLUDEDIRS) $(CCFLAGS) $(LFLAGS) 

%.o: %.c
	$(CC) $(CFLAGS) -c $*.c $(INCLUDEDIRS) -o $*.o

%.o: %.cc
	$(CXX) $(CCFLAGS) -c $*.cc $(INCLUDEDIRS) -o $*.o

clean:
	$(RM) *.o *.exe *.pdf *.tre *.txt *.so *~









