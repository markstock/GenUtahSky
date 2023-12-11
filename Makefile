CC=cc
CFLAGS=-std=c99
LIBS=-lm
ifdef LIBNOVA
  DEFS=-DLIBNOVA
  LIBS+=-lnova
endif

all : genutahsky getsunvec

% : %.c astronomy.c
	$(CC) $(CFLAGS) $(DEFS) -o $@ $^ $(LIBS)

clean :
	rm -f *.o genutahsky getsunvec
