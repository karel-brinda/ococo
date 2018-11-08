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
HTSLIB_VERSION = b6aa0e6

ofiles    = src/main.cpp.o src/misc.cpp.o src/params.cpp.o
hfiles    = $(wildcard src/*.h)

.PHONY: all clean install ococo readme

all: ococo

install: ococo
	install  ococo $(BINDIR)/ococo
	install  $(MANPAGE) $(MANDIR)/$(MANPAGE)

ococo: $(HTSLIB) $(ofiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) $(ofiles) -o $@ -L. $(LIBS) $(HTSLIB)

src/%.cpp.o: src/%.cpp $(hfiles)
	$(CXX) $(CXXFLAGS) $(DFLAGS) -c $< -I $(HTSLIBINCLUDE) -o $@

$(HTSLIB): $(HTSLIBDIR)/Makefile
	$(MAKE) -C $(HTSLIBDIR) lib-static

$(HTSLIBDIR)/Makefile:
	 cd $(HTSLIBDIR) && curl -L https://github.com/samtools/htslib/archive/$(HTSLIB_VERSION).tar.gz | tar xz --strip-components 1

readme:
	f=$$(mktemp);\
	  echo $$f;\
	  sed '/USAGE-BEGIN/q' README.md >> $$f; \
	  printf -- '-->\n```' >> $$f; \
	  man ./ococo.1 | col -b | grep -A999999 SYNO | grep -B99999999 AUTH | ghead -n -1 >> $$f; \
	  printf '```\n<!---\n' >> $$f; \
	  sed -n '/USAGE-END/,$$ p' README.md >> $$f;\
	  cat $$f \
	  | perl -pe 's/^[\s]+$$/\n/g' \
	  | perl -pe 's/[\s]+$$/\n/g' \
	  > README.md
	markdown_py README.md > README.html

clean:
	$(MAKE) -C ext/htslib clean
	rm -f src/*.o
	rm -f ococo

