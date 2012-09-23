# Author Shlomi Lifshits (lifshits_il@yahoo.com)
# Date: September 2012 

CC = gcc
TEMPDIR := $(shell mktemp)

.PHONY: test

all: average_NA.so

test: average_NA.so
	R --vanilla < test.R > $(TEMPDIR)
	cat $(TEMPDIR)

average_NA.so: average_NA.c average_NA.h
	R CMD COMPILE average_NA.c
	R CMD SHLIB average_NA.o

average_NA.c: average_NA.h

clean:
	rm average_NA average_NA.o average_NA.so