#ifndef HL_HUNTER_HH
#define HL_HUNTER_HH

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

class Triple 
{
public: 
  int first_s;
  int second_s;
  int common_length;
  Triple(int f, int s, int c);   
};

std::ostream& operator<<(std::ostream& out, const Triple& T);

class Hunter {
public: 
  Hunter(std::string filename);   // Set up a hunter object from the specified file name (created with translate.cpp)
  ~Hunter() {fp.close();}
  int getWidth() const;           // Get the width of the PH model employed in translate
  int getModel() const;           // Get the actual model employed in translate
  std::list<int> getPositions(std::string seq);  // Get all start positions of the sequence
  std::list<Triple> getTriples(std::string s);   // Gets start positions, combines overlaps.  (Only for 111...1 model.)
  
private:
  static const int int_size = sizeof(int);
  static const int int_size2 = sizeof(int) << 1;
  std::string filename;
  std::ifstream fp;
  int width;
  int model;
  int vec_len;
  int reserve;
  int header_space;
};

inline int translation(const char c) {
  switch(toupper(c)){
  case 'A':{return 0; break;}
  case 'C':{return 1; break;}
  case 'G':{return 2; break;}
  case 'T':{return 3; break;}
  default: {return 4; break;}
  }   
}

inline int translate_seq(const char* seq) {
  int code = 0;
  for (int i=0; seq[i] != '\0'; i++) {
    int t = translation(seq[i]);
    if (t == 4)
      return -1;
    else
      code = (code << 2) | translation(seq[i]);
  }
  return code;
}

#endif

