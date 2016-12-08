OBJECTS = basenode.o node.o gnode.o substmodel.o peel.o userUnif.o RngStream.o io.o reduceSeqs.o
INCLUDEDIRS = $(shell R CMD config --cppflags) $(shell pkg-config --cflags libRmath libR) 
INCLUDEDIRS += -I. -I./models/ehiAIDS
#INCLUDEDIRS += -I$(MKL_INC)
#INCLUDEDIRS += -I. -I./models/ehi3
CCFLAGS = -fopenmp -g -DDEBUG 
CFLAGS = -std=gnu99 -fpic -g -O2 
#LFLAGS = $(shell R CMD config --ldflags) -lgfortran -lquadmath 
#LFLAGS = $(shell R CMD config LAPACK_LIBS)
#LFLAGS += $(shell R CMD config BLAS_LIBS)
LFLAGS = $(shell pkg-config --libs libRmath)
LFLAGS += -lgfortran -lquadmath 
LFLAGS += $(shell pkg-config --libs libR)
#LFLAGS += -lgfortran -lquadmath -lR
GXX = icpc #-fp-model precise #-fp-trap="invalid,divzero,underflow,overflow"
GCC = icc -mkl #-fp-model precise #-fp-trap="invalid,divzero,underflow,overflow"
#GXX = g++
#GCC =  gcc


REXE = R --vanilla

default: genpomp.so

%.so: %.c	
	$(REXE) CMD SHLIB $*.cc -o $*.so

%.exe : ../standalone/%.cc
	$(GXX) -o ../standalone/$*.exe ../standalone/$*.cc $(OBJECTS) $(INCLUDEDIRS) $(CCFLAGS) $(LFLAGS)

% : unit_tests/%.cc
	$(GXX) -mkl -o unit_tests/$* unit_tests/$*.cc $(OBJECTS) $(INCLUDEDIRS) $(CCFLAGS) $(LFLAGS)

%.exe : %.cc $(OBJECTS)
	$(GXX) -o $*.exe $*.cc $(OBJECTS) $(INCLUDEDIRS) $(CCFLAGS) $(LFLAGS)

%.o: %.cc
	$(GXX) $(CCFLAGS) -c $*.cc $(INCLUDEDIRS) -o $*.o

%.o: %.c
	$(GCC) $(CFLAGS) $(INCLUDEDIRS) -c $*.c -o $*.o

clean:
	$(RM) *.o *.exe *.pdf *.tre *.txt *.so 
