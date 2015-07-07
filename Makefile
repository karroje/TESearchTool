OPTS = -g

DIRS = common_code hunter #pHMM_code

all:
	for d in ${DIRS}; do cd $$d; echo $$d; make all; cd ..; done

hunter/hunter.o: common_code/libseq.a force_look
	cd hunter; ${MAKE} hunter.o OPTS=${OPTS}

common_code/libseq.a: force_look
	cd common_code; ${MAKE} OPTS=${OPTS

clean:
	for d in ${DIRS}; do cd $$d; echo $$d; make clean; cd ..; done

force_look:
	true
