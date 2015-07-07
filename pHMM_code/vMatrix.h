#ifndef _V_MATRIX_CLASS
#define _V_MATRIX_CLASS

#include "constants.hh"
#include "matrix.hh"
#include <cmath>

/*
 Specialized Group of Matrices for viterbi searches with a pHMM

 This is the class thats really needed, and is just facilitated by the Matrix<type> class using
 the boost matrix class in my vMatrix implementation is probably the way to go
 */

class vMatrix{
public:

 enum mType{mMatch, mInsert, mDelete};

  vMatrix(int m , int n){
    // this->len = m;
    this->states = n;
    //   Match  = Matrix<double>(len+1,states+1);
    // Insert = Matrix<double>(len+1,states+1);
    //Delete = Matrix<double>(len+1,states+1);
    // End = std::vector<double>(len);
    this->init(m);
  }

  vMatrix(){

  }

  Matrix<double> getMatrix(mType type){
    if (type == mMatch){
      return Match;
    }
    else if (type == mInsert){
      return Insert;
    }
    else if (type == mDelete){
      return Delete;
    }
  }
  
  void init(int l){
    this->len = l;
    Match         = Matrix<double>(len,states);
    MatchPointer  = Matrix<state_type>(len,states);

    Insert        = Matrix<double>(len,states);
    InsertPointer = Matrix<state_type>(len,states);

    Delete        = Matrix<double>(len,states);
    DeletePointer = Matrix<state_type>(len,states);

    End        = vector<double>(len);
    EndPointer = vector<int>(len);

    maxEndPointerIndex = -1;

    for (int j=0; j < len; ++j){
      for (int k = 0; k < states; ++k){
	Match[j][k] = neg_inf;
	Insert[j][k] = neg_inf;
	Delete[j][k] = neg_inf;
      }
      End[j] = neg_inf;
      EndPointer[j] = -1;
    }
  }

  Matrix<double> Match;
  Matrix<state_type> MatchPointer;

  Matrix<double> Insert;
  Matrix<state_type> InsertPointer;

  Matrix<double> Delete;
  Matrix<state_type> DeletePointer;

  vector<double> End;
  vector<int> EndPointer;

  int maxEndPointerIndex;
  
  int size(){
    return len;
  }

  int numStates(){
    return states;
  }

private:
  int len;
  int states;

protected:

};
#endif
