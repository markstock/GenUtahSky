CC=cc
CFLAGS=-std=c99 -Wall
LIBS=-lm
ifdef LIBNOVA
  DEFS=-DLIBNOVA
  LIBS+=-lnova
endif

all : genutahsky getsunvec

% : %.c astronomy.c
	$(CC) $(CFLAGS) $(DEFS) -o $@ $^ $(LIBS)

install :
	lastpath=$(echo "${RAYPATH}" | tr ':' '\n' | tail -n 1)
	[[ -z "${lastpath}" ]] && echo "Error: RAYPATH is not defined; set it and re-run 'make install'"
	[[ -n "${lastpath}" ]] && cp *.cal ${lastpath} && cp TychoSkymapII*hdr ${lastpath} && cp stardome.rad ${lastpath}

clean :
	rm -f *.o genutahsky getsunvec
