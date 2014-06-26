include makefile.inc

SUBDIRS=src

all: recursive
	for i in $(SUBDIRS); do cd "$$i" && $(MAKE) $@ && cd ..; done

recursive:
	true

clean:
	for i in $(SUBDIRS); do \
		cd "$$i" && $(MAKE) clean && cd ..; \
	done
