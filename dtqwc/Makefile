CXX=g++
CXXFLAGS=-std=gnu++11 -O3 -Wl,--no-as-needed -m64 -march=native -mtune=native -fopenmp -I/usr/local/include -I/usr/local/include/eigen3

DEPS=

ODIR=obj

OBJnewdtq=newdtq.o
newdtqOBJ=$(patsubst %,$(ODIR)/%,$(OBJnewdtq))

OBJnewdtqinf=newdtqinf.o
newdtqinfOBJ=$(patsubst %,$(ODIR)/%,$(OBJnewdtqinf))

OBJdtqceres=dtqceres.o
dtqceresOBJ=$(patsubst %,$(ODIR)/%,$(OBJdtqceres))

LDFLAGS=-L/opt/intel/mkl/lib/intel64
LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp -larmadillo
NLOPT=-L/usr/local/lib -lnlopt
CERES=-L/usr/local/lib64 -lceres -L/usr/local/lib -lglog

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

newdtq: $(newdtqOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

newdtqinf: $(newdtqinfOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(NLOPT)

dtqceres: $(dtqceresOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS) $(CERES)

.PHONY: clean

clean:
	rm -f newdtq newdtqinf $(ODIR)/*.o *~ core 

