OBJECTS = basenode.o node.o gnode.o substmodel.o tree.o peel.o userUnif.o RngStream.o io.o reduceSeqs.o
INCLUDEDIRS = $(shell R CMD config --cppflags) $(shell pkg-config --cflags libRmath libR) 
CCFLAGS = -fopenmp -g -D DEBUG
CFLAGS = -std=gnu99 -fpi -g #-O2
#LFLAGS = $(shell R CMD config --ldflags) -lgfortran -lquadmath 
LFLAGS = $(shell R CMD config LAPACK_LIBS)
LFLAGS += $(shell R CMD config BLAS_LIBS)
LFLAGS += $(shell pkg-config --libs libRmath)
LFLAGS += -lgfortran -lquadmath 
LFLAGS += $(shell pkg-config --libs libR)
#LFLAGS += -lgfortran -lquadmath -lR

%.exe: %.cc
	cd ../src; make $(OBJECTS)
	cd ../src; make $*.exe

#%.exe : %.cc ../$(OBJECTS)
#	g++ -o $*.exe $*.cc $(OBJECTS) $(INCLUDEDIRS) $(CCFLAGS) $(LFLAGS)

%.o: %.cc
	g++ $(CCFLAGS) -c $*.cc $(INCLUDEDIRS) -o $*.o

%.o: %.c
	gcc $(CFLAGS) $(INCLUDEDIRS) -c $*.c -o $*.o

clean:
	$(RM) *.o *.exe *.pdf *.tre *.txt *~
	cd ../src; make clean
