all : genutahsky getsunvec

genutahsky : genutahsky.c
	cc -o genutahsky genutahsky.c -lm -lnova

getsunvec : getsunvec.c
	cc -o getsunvec getsunvec.c -lm -lnova
