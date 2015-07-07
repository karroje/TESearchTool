#ifndef HL_HUNTER_HH
#define HL_HUNTER_HH

#include <list>
#include <string>
#include <iostream>
#include "HashTable.h"
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
  Hunter(std::string fn,int width);
  void buildHash(); 
  std::list<std::string> readFasta(string seqfile);
  int getWidth()const;
  int getModel()const;
  std::list<Triple> getTriples(std::string s);  
  
private:
  static const int int_size = sizeof(int);
  std::string filename;
  int width;        
  std::list<HashTable*>::const_iterator hash_it;//Keeps track of the fragment in use
  std::list<HashTable*> hashes;//Hashtables for each fragment
};

#endif

