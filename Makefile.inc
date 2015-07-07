CC = g++
DEBUG = -g
OPTS = ${DEBUG} -Wall -Wno-sign-compare

COMMON_CODE = -I../ -I../common_code
SEQL = -L../common_code -lseq
SEQUENCEL = ${SEQL} -L../common_code -lsequence
INC = -I/opt/local/include/ 

BOOST_INC = /opt/local/include/
BOOST_LIB = /opt/local/lib/
BOOST = -I${BOOST_INC} -L${BOOST_LIB}

HUNTER = ../hunter/hunter.o

