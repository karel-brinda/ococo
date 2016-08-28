.PHONY: all clean

all: misc.o main.o params.o

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -I ../ext/htslib/

clean:
	rm -f *.o
