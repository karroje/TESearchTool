include ../Makefile.inc

all: translate TestHunter hunter.o HashTable.o Tester.o

translate: translate.cpp hunter.o ../common_code/libseq.a
	${CC} -o translate translate.cpp hunter.o ${SEQL} ${OPTS} ${INC} ${COMMON_CODE}

TestHunter:tester.cpp hunter.o HashTable.o
	${CC} -o TestHunter hunter.cpp hunter.o HashTable.cpp HashTable.o${SEQL} ${OPTS} ${INC} ${COMMON_CODE}

Tester.o: Tester.cpp Tester.h 
	${CC} -c Tester.cpp ${SEQL} ${OPTS} ${INC} ${COMMON_CODE} 

hunter.o: hunter.cpp hunter.hh HashTable.h
	${CC} -c hunter.cpp ${SEQL} ${OPTS} ${INC} ${COMMON_CODE} 

HashTable.o: HashTable.h
	${CC} -c HashTable.cpp ${SEQL} ${OPTS} ${INC} ${COMMON_CODE} 

../common_code/libseq.a:
	cd ../common_code/libseq.a; make libseq.a


clean:
	rm *.o translate