// endProbFunc.hh
#include "constants.hh"
#include <iostream>
#include <string>
#include <cstdlib>

#ifndef END_PROB_FUNC_HH
#define END_PROB_FUNC_HH


class endProbFunc {
public:
  endProbFunc(int num_states, int left_anchor, int right_anchor) {
    this->num_states = num_states;
    this->left_anchor = left_anchor;  
    this->right_anchor = right_anchor;
  }

  virtual double getProb(int state) = 0;

protected:
  int num_states;
  int left_anchor;
  int right_anchor;
};

class fixedEndpointDist : public endProbFunc {
public:
  fixedEndpointDist(int num_states, int left_anchor, int right_anchor) : endProbFunc(num_states, left_anchor, right_anchor) {}
  virtual double getProb(int state) {
    if (state == 0 || state == num_states-1)
      return 0;
    return neg_inf;
  }
};

class uniformEndpointDist : public endProbFunc {
public:

  // Formulats reflect the pmf for Z = min(X,Y) and Z=max(X,Y), X,Y uniform.
  // Each state number adjusted to count from 1 instead of 0.
  uniformEndpointDist(int num_states, int left_anchor, int right_anchor) : endProbFunc(num_states, left_anchor, right_anchor) {}
  virtual double getProb(int state) {
    if (state <= left_anchor) 
      return log((double)2*(num_states - state - 1)) - log((double)2*num_states*(left_anchor+1) / ((left_anchor+1)*(left_anchor+1)));

    else if (state >= right_anchor) 
      return log((double)2*state+1) - log((double)(num_states - right_anchor)) - log((double)(num_states + right_anchor));

    std::cerr << "UniformEndpointDist:Bad state parameters\n";
    exit(1);
    return 0;
  }
};
      

class endProbFuncFactory {
public:
  endProbFuncFactory();
  endProbFuncFactory(std::string probType);
  void setType(std::string probType);
  endProbFunc* createEndProbFunc(int num_states, int left_anchor, int right_anchor) const;
private:
  int probType;
};

#endif // END_PROB_FUNC_HH
