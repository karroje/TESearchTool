#include "alignment.hh"
#include <string>
#include "seq.h"

using namespace std;

pair<string, string> Alignment::reconstruct(SEQ* modern, string ancestor) {
  uchar *m = modern->seq + m_start;
  string::iterator a = ancestor.begin() + a_start;
  
  int s = size();
  string modern_aligned(s-2, 0);
  string ancestor_aligned(s-2, 0);

  string::iterator mi = modern_aligned.begin();
  string::iterator ai = ancestor_aligned.begin();

  for (stateVec_iterator i=stateType.begin(); i != stateType.end(); i++) {
    switch (*i) {
    case mat: *mi = *m; m++; *ai = *a; a++; mi++; ai++; break;
    case ins: *mi = *m; m++; *ai = '-'; mi++; ai++; break;
    case del: *mi = '-'; *ai = *a; a++; mi++; ai++; break;
    case end: break;					 
    }
  }

  return make_pair(modern_aligned, ancestor_aligned);
}

void Alignment::print_alignment(ostream& out, SEQ* modern, string ancestor, int col_width) {
  pair<string,string> p = reconstruct(modern, ancestor);
  int s = p.first.size();

  for (int i=0; i <= s/col_width; i++) {
    out << p.first.substr(col_width*i, col_width) << "\n";
    for (int j=col_width*i; j < col_width*(i+1) && j < (int)p.first.size(); j++) {
      switch (stateType[j+1]) {
      case mat: out << ((p.first[j]==p.second[j]) ? "|" : "."); break;
      case ins: out << " "; break;
      case del: out << " "; break;
      case end: break;
      }
    }
    out << "\n";

    out << p.second.substr(col_width*i, col_width) << "\n\n";
  }
}
  

Alignment combine_alignment(Alignment& L, Alignment& R, int anchor_length) {
  if (L.score==neg_inf || R.score==neg_inf)
    return Alignment();

  Alignment A(L.size() + R.size() + anchor_length - 2);
  stateVec_iterator p = copy(L.stateType.begin(), L.stateType.end(), A.stateType.begin());
  for (int i=1; i < anchor_length-1; i++)
    *p++ = mat;
  copy(R.stateType.begin(), R.stateType.end(), p);

  A.score = L.score + R.score;
  A.a_start = L.a_start;
  A.m_start = L.m_start;

  return A;
}
  
