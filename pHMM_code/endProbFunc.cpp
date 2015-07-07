// endProbFun.cpp
#include "endProbFunc.hh"
#include <string>


using namespace std;



endProbFuncFactory::endProbFuncFactory() {
  this->probType = -1;
}

endProbFuncFactory::endProbFuncFactory(string probType) {
  setType(probType);
}

void endProbFuncFactory::setType(string probType) {
  if (probType == "fixed")
    this->probType = 0;
  if (probType == "uniform") 
    this->probType = 1;
  else {
    cerr << "endProbFuncFactory::setType given bad parameteter: " << probType << endl;
    exit(1);
  }
}

endProbFunc* endProbFuncFactory::createEndProbFunc(int size, int left_anchor, int right_anchor) const {
 switch (probType) {
 case 0: return new fixedEndpointDist(size, left_anchor, right_anchor);
 case 1: return new uniformEndpointDist(size, left_anchor, right_anchor);
 default: cerr << "endProbFuncFactory object used without proper initiliziation.";
 }
 return NULL;
}
					
