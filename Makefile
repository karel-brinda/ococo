CXX      ?= g++
CXXFLAGS  = -std=c++11 -Wall -Wextra -Wno-missing-field-initializers -Wshadow -g -O2
LIBS      = -lm -lz -lpthread

PREFIX    = $(DESTDIR)/usr/local
BINDIR    = $(PREFIX)/bin
MANDIR    = $(PREFIX)/share/man/man1
MANPAGE   = ococo.1

HTSLIBDIR = ext/htslib
HTSLIB    = $(HTSLIBDIR)/libhts.a

export CXX
export CXXFLAGS

export HTSLIBDIR
export HTSLIB

.PHONY: all clean install ococo

all: ococo

install: ococo
	install  ococo $(BINDIR)/ococo
	install  $(MANPAGE) $(MANDIR)/$(MANPAGE)

ococo: $(HTSLIB)
	$(MAKE) -C ./src

	$(CXX) $(CXXFLAGS) $(DFLAGS) ./src/*.o -o $@ -L. $(LIBS) $(HTSLIB)

$(HTSLIB):
	$(MAKE) -C $(HTSLIBDIR) lib-static

clean:
	$(MAKE) -C ext/htslib clean
	$(MAKE) -C src clean
	rm ococo
