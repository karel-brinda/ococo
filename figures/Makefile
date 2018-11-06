.PHONY: all help clean cleanall

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

all: Figure_1.png

Figure_1.png: Figure_1.pdf Makefile
	convert -flatten -density 300 "$<" -resize 50% "$@"

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean

cleanall: clean ## Clean all


