// dna.h
#include <string>
#pragma once

enum nucleotide {_A, _C, _G, _T, _R, _Y, _N};

const int num_base_codes = 7;

inline nucleotide getNucleotide(char x) {
  switch (x) {
  case 'A':
  case 'a': return _A;
  case 'C':
  case 'c': return _C;
  case 'G':
  case 'g': return _G;
  case 'T':
  case 't': return _T;
  case 'R': 
  case 'r': return _R;
  case 'Y':
  case 'y': return _Y;
  }
  return _N;
}

inline char convertNucleotide(nucleotide n) {
  switch (n) {
  case _A: return 'A';
  case _C: return 'C';
  case _G: return 'G';
  case _T: return 'T';
  case _R: return 'R';
  case _Y: return 'Y';
  case _N: return 'N';
  }

  return 'N';
}

inline char converNucleotide(int n) {
  return convertNucleotide((nucleotide)n);
}

