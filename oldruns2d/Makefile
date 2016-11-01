CXX = g++
CXXFLAGS = -I/usr/include/R -I/home/hbhat/R/x86_64-redhat-linux-gnu-library/3.1/Rcpp/include -I/home/hbhat/R/x86_64-redhat-linux-gnu-library/3.2/RcppArmadillo/include -O3 -m64 -march=native -mtune=native -fopenmp
LDFLAGS =
LIBS = -lgsl -lgslcblas -lm -lgomp

TARGETS = 

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


