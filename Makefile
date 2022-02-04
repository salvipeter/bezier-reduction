all: reduction-test

CFLAGS=-std=c11 -Wall -pedantic -O3
LIBS=-lgmp

reduction-test: reduction-test.o reduction.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
