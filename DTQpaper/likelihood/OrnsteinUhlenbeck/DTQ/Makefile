CXX=g++
CXXFLAGS=-std=c++1y -O3 -Wl,--no-as-needed -m64 -march=native -mtune=native -fopenmp -I/usr/local/include 

DEPS=

ODIR=obj

OBJdtqou10sweep=dtqou10sweep.o
dtqou10sweepOBJ=$(patsubst %,$(ODIR)/%,$(OBJdtqou10sweep))

OBJdtqou20sweep=dtqou20sweep.o
dtqou20sweepOBJ=$(patsubst %,$(ODIR)/%,$(OBJdtqou20sweep))

LDFLAGS=-L/opt/intel/mkl/lib/intel64
LIBS=-lmkl_intel_lp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lgomp -larmadillo

$(ODIR)/%.o: %.cpp $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

dtqou10sweep: $(dtqou10sweepOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

dtqou20sweep: $(dtqou20sweepOBJ)
	$(CXX) -o $@ $^ $(CXXFLAGS) $(LDFLAGS) $(LIBS)

.PHONY: clean

clean:
	rm -f dtqou10sweep dtqou20sweep $(ODIR)/*.o *~ core 

