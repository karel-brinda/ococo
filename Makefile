CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -g -O2
LIBS      = -lm -lz -lpthread

PREFIX    = $(DESTDIR)/usr/local
BINDIR    = $(PREFIX)/bin
MANDIR    = $(PREFIX)/share/man/man1
MANPAGE   = ococo.1

HTSLIBDIR = ext/htslib
HTSLIB    = $(HTSLIBDIR)/libhts.a
HTSLIBINCLUDE = $(HTSLIBDIR)

ofiles    = src/main.o src/misc.o src/params.o
hfiles    = $(wildcard src/*.h)

.PHONY: all clean install ococo

all: ococo

install: ococo
	install  ococo $(BINDIR)/ococo
	install  $(MANPAGE) $(MANDIR)/$(MANPAGE)

ococo: $(HTSLIB) $(ofiles) 
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(ofiles) -o $@ -L. $(LIBS) $(HTSLIB)

src/%.o: src/%.cpp src/%.h $(hfiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -I $(HTSLIBINCLUDE) -o $@

$(HTSLIB):
	$(MAKE) -C $(HTSLIBDIR) lib-static

clean:
	$(MAKE) -C ext/htslib clean
	rm -f src/*.o
	rm -f ococo

