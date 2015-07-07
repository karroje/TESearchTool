/*
	@author Sean Wilkerson  <wilkersm@muohio.edu>
*/
#include "alignment.hh"
#include "matrix.hh"
#include "vMatrix.h"
#include "constants.hh"
#include "endProbFunc.hh"
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <limits>
#include <boost/tuple/tuple.hpp>
#include "seq.h"
#include "dna.h"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/math/distributions/binomial.hpp>

#include <map>
#ifndef pHMM_h
#define pHMM_h

using namespace std;


enum char_map{A,C,G,T,R,Y,N};
const int LEFT = -1;
const int RIGHT = 1;

class pHMM{
public:
  //**************************************
  // Constructors / destructors
  pHMM () {};
  pHMM (int m) {resize(m);}
  pHMM (string name, int m) {this->name = name; resize(m);}
  pHMM (char* hmm_file);
  ~pHMM();

  //**************************************
  // Interface methods
  // transProb and emissionProb turned out to be unncessary, but didn't hurt anything
  bool printToFile(string filename, bool logValues = true);
  double length() {return M;}

  // Accessors: get the transition (log) probabilitites
  string getName() {return this->name;}
  string getAncestorSequence() {return this->ANCESTOR;}
  double getAncestralBG(char c) {return ancestral_bg[getNucleotide(c)];}
  int getModelLength() {return this->M;}

  inline double transProb(state_type source, state_type dest) {return t_P[source][dest];}
  inline double transProb(int source, int dest) {return transProb( (state_type) source, (state_type) dest);}
  inline double emissionProb(int col, int obs, state_type stype) {
    switch(stype) {
    case mat: return m_eP[col][obs]; 
    case ins: return i_eP[col][obs];
    case del:
    case end: cerr << "Bad state type passed to pHMM::emissionProb\n"; exit(1);
    }
    return 0;
  }
  inline double emissionProb(int col, char obs, state_type stype) {
    return emissionProb(col, getNucleotide(obs), stype);
  }
  inline double emissionProb(int col, nucleotide obs, state_type stype) {
    return emissionProb(col, (int)obs, stype);
  }

  // Modifiers: set the transition (log) probabilities
  void resize(int m);
  void setName(const string& name) {this->name = name;}
  void setAncestor(const string& ancestor);  // Also sets ancestral_bg

  void setTransProb(state_type source, state_type dest, double value) {t_P[source][dest] = value;}
  void setTransProb(int source, int dest, double value) {setTransProb( (state_type) source, (state_type) dest, value);}

  void setEmissionProb(int col, int obs, state_type stype, double value) {
    switch(stype) {
    case mat: m_eP[col][obs] = value; return;
    case ins: i_eP[col][obs] = value; return;
    case del:
    case end: cerr << "Bad state type passed to pHMM::emissionProb\n"; exit(1);
    }
    return;
  }
  void setEmissionProb(int col, char obs, state_type stype, double value) {
    setEmissionProb(col, getNucleotide(obs), stype, value);
  }
  void setEmissionProb(int col, nucleotide obs, state_type stype, double value) {
    setEmissionProb(col, (int)obs, stype, value);
  }


  //**************************************
  // Automated construction methods
  // void buildFromAncestor(vector<string> &filedump,string &name);


  //**************************************
  // HMM Algorithms
  // bDviterbi: if "use_insertion" == false, then background distribution will be used for insertion distribution.
  //            Otherwise, the HMM emission probability for insertion is used
  // bDforward:  
  Alignment bDviterbi(SEQ* obs, int start_left, int start_right, int rel_left, int rel_right, double* normalized_bg_freq, endProbFunc* EPFunc, bool use_insertion = false); //bi-directional modified viterbi
  boost::tuple<int, int, double> bDforward(SEQ* obs, int start_left, int start_right, int rel_left, int rel_right, double* normalized_bg_freq, endProbFunc* EPFunc, bool use_insertion);

  //************************************
  //Generating Sequences
  pair<string,string> genSequence();    // KARRO function	
  string genBackground(int bg_length);


private:
  //**************************************
  // Attributes
  string name;
  string ANCESTOR;
  vector<double> ancestral_bg; // estimated background distribution of ancestor 

  Matrix<double> i_eP; // insert emission probs
  Matrix<double> m_eP; // match emission probs
  Matrix<double> t_P; // transition probs
  vector<double> end_P; // ending probs 

  int M; // model length

  // For use in HMM algorithms
  vMatrix left;
  vMatrix right;

  // Characteristics of the alphabet
  void build(Matrix<char>&);

  //**************************************
  // helper methods 
  void initAlphabet(string alphabet = "ACGTRYN");
    
  //other misc helper stuff
  bool inBounds(int i, int stop, int dir);
  state_type getStateType(char base, char ancestor); //dependent upon .mal file formatting


  //helper functions for bidirectional search
  void bDv_fill(SEQ* obs, int start, int end,  vMatrix& dpm, int q_coord, int dir, int o2r, double* normalized_bg_freq, bool use_insertion, endProbFunc* EPFunc);
  Alignment bDv_trace(vMatrix& dpm, int start, int qstate, int dir);

  void bDf_fill(SEQ* obs, int start, int end, vMatrix& dpm, int q_coord, int dir, int o2r, double* normalized_bg_freq, bool use_insertion, endProbFunc* EPFunc);
};


//*******************************
// Other useful functions and contants


#endif

