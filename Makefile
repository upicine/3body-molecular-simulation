#
# A template for the 2016 MPI lab at the University of Warsaw.
# Copyright (C) 2016, Konrad Iwanicki.
#
CC          := mpic++
# Available flags:
# -DUSE_RANDOM_GRAPH=1   --- generates a random graph
# -DUSE_RANDOM_SEED=123  --- uses a given seed to generate a random graph
CFLAGS      := -O3 -c
LFLAGS      := -O3
ALL         := \
	main.exe


all : $(ALL)

main.exe : main.o particle-parser.o utils.o embedded-algorithm.o
	$(CC) $(LFLAGS) -o $@ $^

main.o : main.cpp utils.h particle-parser.h
	$(CC) $(CFLAGS) $<

utils.o : utils.cpp utils.h
	$(CC) $(CFLAGS) $<

particle-parser.o : particle-parser.cpp particle-parser.h utils.h
	$(CC) $(CFLAGS) $<

embedded-algorithm.o : embedded-algorithm.cpp embedded-algorithm.h utils.h
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o *.out *.err $(ALL)
