CC          := mpic++
CFLAGS      := -O3 -Wall -c
LFLAGS      := -O3 -Wall
ALL         := \
	main.exe


all : $(ALL)

main.exe : main.o particle-parser.o utils.o embedded-algorithm.o verlet-integration.o particle-buffer.o
	$(CC) $(LFLAGS) -o $@ $^

main.o : main.cpp utils.h particle-parser.h
	$(CC) $(CFLAGS) $<

utils.o : utils.cpp utils.h
	$(CC) $(CFLAGS) $<

particle-parser.o : particle-parser.cpp particle-parser.h utils.h
	$(CC) $(CFLAGS) $<

embedded-algorithm.o : embedded-algorithm.cpp embedded-algorithm.h verlet-integration.h utils.h particle-buffer.h
	$(CC) $(CFLAGS) $<

verlet-integration.o : verlet-integration.cpp utils.h verlet-integration.h embedded-algorithm.h particle-buffer.h
	$(CC) $(CFLAGS) $<

particle-buffer.o : particle-buffer.cpp particle-buffer.h utils.h
	$(CC) $(CFLAGS) $<

clean :
	rm -f *.o *.out *.err $(ALL)
