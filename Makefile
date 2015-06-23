all : genutahsky

genutahsky : genutahsky.c
	cc -o genutahsky genutahsky.c -lm -lnova
