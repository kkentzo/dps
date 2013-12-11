CC=gcc
CFLAGS=-Wall -g # -pg
LDFLAGS=-lgsl -lgslcblas -lhdf5 -lhdf5_hl `pkg-config --cflags --libs glib-2.0`


dps: *.c lib/clib.c lib/sampling.c lib/distribution.c lib/nvar.c lib/qhdf.c lib/bitvector.c
