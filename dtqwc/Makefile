CXX=g++
CXXFLAGS=-std=gnu++11 -O3 -Wl,--no-as-needed -m64 -march=native -mtune=native -fopenmp

DEPS=

ODIR=obj

OBJnewdtq=newdtq.o
newdtqOBJ=$(patsubst %,$(ODIR)/%,$(OBJnewdtq))

OBJnewdtqinf=newdtqinf.o
newdtqinfOBJ=$(patsubst %,$(ODIR)/%,$(OBJnewdtqinf))

LDFLAGS=-L/opt/intel/mkl/lib/intel64
LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp -larmadillo
NLOPT=-L/usr/local/lib -lnlopt

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

newdtq: $(newdtqOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

newdtqinf: $(newdtqinfOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(NLOPT)

.PHONY: clean

clean:
	rm -f newdtq newdtqinf $(ODIR)/*.o *~ core 

