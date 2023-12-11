CC=cc
CFLAGS=-std=c99
DEFS=-DLIBNOVA -DASTRO
LIBS=-lm -lnova

all : genutahsky getsunvec

% : %.c astronomy.c
	$(CC) $(CFLAGS) $(DEFS) -o $@ $^ $(LIBS)

clean :
	rm -f *.o genutahsky getsunvec
