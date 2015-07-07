include ../Makefile.inc

all: libseq.a libsequence.a

libseq.a: seq.c nib.c util.c charvec.c
	gcc -c seq.c nib.c util.c charvec.c ${OPTS}
	ar rc libseq.a seq.o nib.o util.o charvec.o 
	ranlib libseq.a
	rm seq.o nib.o util.o charvec.o 

libsequence.a: sequence.cpp
	g++ -c sequence.cpp ${OPTS}
	ar rc libsequence.a sequence.o
	ranlib libsequence.a
	rm sequence.cpp

clean:
	rm libseq.a libsequence.a