CXX      ?= g++
CXXFLAGS = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -Wshadow -g -O2
LIBS     = -lm -lz -lpthread

export CXX
export CXXFLAGS

.PHONY: all clean

all: ococo

ococo:
	$(MAKE) -C ./ext/htslib lib-static
	$(MAKE) -C ./src

	$(CXX) $(CXXFLAGS) $(DFLAGS) ./src/*.o -o $@ -L. $(LIBS) ./ext/htslib/libhts.a

ext/htslib/libhts.a:
	$(MAKE) -C ext/htslib lib-static

clean:
	$(MAKE) -C ext/htslib clean
	$(MAKE) -C src clean
	rm ococo
