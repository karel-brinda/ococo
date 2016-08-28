.PHONY: all clean

all: misc.o main.o params.o

HEADERS=types.h stats.h version.h

main.o: main.cpp consensus.h caller.h version.h $(HEADERS)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -I ../$(HTSLIBINCLUDE)


%.o: %.cpp %.h $(HEADERS)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -I ../$(HTSLIBINCLUDE)

clean:
	rm -f *.o
