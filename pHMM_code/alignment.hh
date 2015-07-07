#include "seq.hh"
#include "vMatrix.h" // Needed to get the neg_inf def
#include <fstream>

#ifndef JEK_ALIGNMENT_HH
#define JEK_ALIGNMENT_HH

typedef vector<state_type> stateVec;
typedef vector<state_type>::iterator stateVec_iterator;

class Alignment {
public:
  Alignment() {
    this->score = neg_inf;
  }

  Alignment(int size) {
    this->score = neg_inf;
    this->stateType = stateVec(size);
  }

  Alignment(double score, int ancestor_start, int modern_start, stateVec stateType = stateVec()) {
    this->score = score;
    this->a_start = ancestor_start;
    this->m_start = modern_start;
    this->stateType = stateType;
  }


  int size() const {return stateType.size();}
  std::pair<int,int> modern_ends() {return std::make_pair(m_start, m_start + size() - 2);}
  std::pair<int,int> ancestor_ends() {return std::make_pair(a_start, a_start + size() -2);}

  bool isEmptyPath() {return this->score == neg_inf;}

  pair<string,string> reconstruct(SEQ* modern, string ancestor);
  void print_alignment(ostream& out, SEQ* moden, string ancestor, int col_width = 80);

  double score;
  int a_start;             // start coordinate of the Alignment on the ancestor (= HMM column number)
  int m_start;        // start coordinate of the Alignment on the genome
  stateVec stateType;   // the Alignment state-type sequence
};


Alignment combine_alignment(Alignment& L, Alignment& R, int anchor_size);

#endif
