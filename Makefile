CC              = g++ -std=c++11
CC_OBJ_FLAGS    = -c -fPIC
CC_Shared_FLAGS = -shared -Wl,-soname,libJPsiGen.so

ROOTCONFIG      = /cvmfs/oasis.opensciencegrid.org/jlab/scicomp/sw/el9/root/6.30.02-gcc11.4.0/bin/root-config
ROOT_CFLAGS     = $(shell $(ROOTCONFIG) --cflags)
ROOT_LIBS       = $(shell $(ROOTCONFIG) --libs)

libJPsiGEN      = libJPsiGen

all: JPsiGen.cc TTCSKine.o KinFuncs.o
	@mkdir -p lib ; rm -f lib/*.so
	$(CC) $(CC_Shared_FLAGS) -o lib/${libJPsiGEN}.so.1.0.1 TTCSKine.o KinFuncs.o
	cd lib && \
	ln -sf ${libJPsiGEN}.so.1.0.1 ${libJPsiGEN}.so.1 && \
	ln -sf ${libJPsiGEN}.so.1.0.1 ${libJPsiGEN}.so
	$(CC) -o JPsiGen.exe JPsiGen.cc -I ./include -L./lib -lJPsiGen $(ROOT_CFLAGS) $(ROOT_LIBS)

TTCSKine.o: src/TTCSKine.cc include/TTCSKine.h
	$(CC) $(CC_OBJ_FLAGS) src/TTCSKine.cc -o $@ $(ROOT_CFLAGS) -I ./include

KinFuncs.o: src/KinFunctions.cc include/KinFunctions.h
	$(CC) $(CC_OBJ_FLAGS) src/KinFunctions.cc -o $@ $(ROOT_CFLAGS) -I ./include

clean:
	rm -f JPsiGen.exe *.o lib/*.so.*

