CXX=g++
CXXFLAGS=-std=gnu++11 -O3 -Wl,--no-as-needed -m64 -march=native -mtune=native -fopenmp

DEPS=

ODIR=obj

_newdtqOBJ=newdtq.o
newdtqOBJ=$(patsubst %,$(ODIR)/%,$(_newdtqOBJ))

_newdtqinfOBJ=newdtqinf.o
newdtqinfOBJ=$(patsubst %,$(ODIR)/%,$(_newdtqinfOBJ))

LDFLAGS=-L/opt/intel/mkl/lib/intel64
LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp -larmadillo
NLOPT=-L/usr/local/lib -lnlopt

$(ODIR)/%.o: %.c $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

newdtq: $(newdtqOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

newdtqinf: $(newdtqinfOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(NLOPT)

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ core 

