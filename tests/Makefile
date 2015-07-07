include ../Makefile.inc

DIR = test1
SIM_PARAMS = sim_params.txt
WORD_SIZE = 3

new_test: reset_test run_test 

run_test: ../pHMM_code/pHMM_TE ${DIR}
	../pHMM_code/pHMM_TE -p ${DIR}/test.hmm -q ${DIR}/test.fa -h ${DIR}/test.hash -e uniform -B 0.1x0.4x0.4x0.1

reset_test: clear_test ${DIR}

${DIR}: ${DIR}/test.fa ${DIR}/test.hash ${DIR}/test.hmm

${DIR}/test.hmm: ../pHMM_code/build_phmm ${DIR}/test.fa 
	../pHMM_code/build_phmm ${DIR}/test.mal

${DIR}/test.hash: ../hunter/translate 
	../hunter/translate ${DIR}/test.fa ${DIR}/test.hash ${WORD_SIZE}

${DIR}/test.fa: ../simulator/basic_simulator.pl ${DIR}/${SIM_PARAMS}
	python ../simulator/basic_simulator.py -c ${DIR}/${SIM_PARAMS} ${DIR}/test.mal ${DIR}/test.fa

clear_test:
	rm -f ${DIR}/test*

clear:
	rm -f ${DIR}/test*