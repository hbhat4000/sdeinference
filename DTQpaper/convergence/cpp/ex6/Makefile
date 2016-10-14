CXX = g++
CXXFLAGS = -O3 -m64 -march=native -mtune=native -fopenmp
LDFLAGS = 
LIBS = -lgsl -lgslcblas -lm -lgomp

TARGETS = iterate
objects = iterate.o

.SUFFIXES: .o .cpp

.coo.o:
	$(CXX) $(CXXFLAGS) -c $*.cpp

%.mod : %.o
	@if [! -f $@ ]; then \
	  rm $< \
	  $(MAKE) $< \
	fi 

$(TARGETS):
	$(CXX) $(LDFLAGS) $(objects) -o $@ $(LIBS)

all:
	$(TARGETS)

iterate:	$(objects)

clean:
	-rm -f *.o $(TARGETS)


