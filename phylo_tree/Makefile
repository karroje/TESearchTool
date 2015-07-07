include ../Makefile.inc

all: translate hunter.o chunks

translate: translate.cpp hunter.o ../common_code/libseq.a
	${CC} -o translate translate.cpp hunter.o ${SEQL} ${OPTS} ${INC} ${COMMON_CODE}


hunter.o: hunter.cpp hunter.hh
	${CC} -c hunter.cpp ${SEQL} ${OPTS} ${INC} ${COMMON_CODE} 

chunks: chunks.cpp hunter.hh
	${CC} -o chunks chunks.cpp ${HUNTER} ${SEQUENCEL} ${OPTS} ${INC} ${COMMON_CODE}

../common_code/libseq.a:
	cd ../common_code/libseq.a; make libseq.a


clean:
	rm hunter.o translate chunks
