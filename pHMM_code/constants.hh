// constants.hh
// Used to define constants for easy reference

#include <math.h>
#include <string>
#pragma once


const double neg_inf = log((double)0);
enum state_type{mat,ins,del,end};

#define INS_CHAR '-' 

#define _TRANS 3

#define MAT_MAT 0
#define MAT_INS 1
#define MAT_DEL 2
#define DEL_DEL 3
#define DEL_MAT 4
#define INS_MAT 5
#define INS_INS 6
#define MAT_END 7

#define STATES 8 // insert,delete,match,end

