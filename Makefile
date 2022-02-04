all: reduction-test

CFLAGS=-std=c11 -Wall -pedantic
LIBS=-lgmp

reduction-test: reduction-test.o reduction.o
	$(CC) -o $@ $^ $(CFLAGS) $(LIBS)
