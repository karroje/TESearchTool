include ../Makefile.inc

all: pHMM_TE build_phmm

pHMM_TE: pHMM_TE.cpp alignment.o endProbFunc.o pHMM.o ${HUNTER}
	${CC} -o pHMM_TE pHMM_TE.cpp alignment.o pHMM.o endProbFunc.o ${HUNTER} ${OPTS} ${COMMON_CODE} ${SEQL} ${INC} ${BOOST}

build_phmm: build_phmm.cpp pHMM.o alignment.o endProbFunc.o
	${CC} -o build_phmm build_phmm.cpp pHMM.o alignment.o endProbFunc.o ${OPTS} ${COMMON_CODE} ${SEQ} ${BOOST}

alignment.o: alignment.cpp alignment.hh
	${CC} -c alignment.cpp ${OPTS} ${COMMON_CODE} ${SEQL}

endProbFunc.o: endProbFunc.cpp endProbFunc.hh
	${CC} -c endProbFunc.cpp ${OPTS} ${COMMON_CODE}

pHMM.o: pHMM.cpp pHMM.hh constants.hh endProbFunc.hh alignment.hh vMatrix.h
	${CC} -c pHMM.cpp ${OPTS} ${INC} ${COMMON_CODE} ${BOOST}

clean:
	rm *.o build_phmm pHMM_TE

