CC=gcc
CFLAGS=-Wall -g # -pg
LDFLAGS=-lgsl -lgslcblas -lhdf5 -lhdf5_hl `pkg-config --cflags --libs glib-2.0`

dps: dps.c params.c logger.c cell.c plasmid.c population.c pool.c ../utils/clib.c ../utils/sampling.c ../utils/distribution.c ../utils/nvar.c ../utils/qhdf.c ../utils/bitvector.c
