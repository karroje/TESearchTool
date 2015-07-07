//#include "pHMM.hh"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctype.h>
#include <assert.h>

using namespace std;

#define DEL_CHAR '-'

enum column_offsets {chr_location, chr_start, chr_stop, ancestor_start, ancestor_stop, modern_sequence};
enum state_type {match, insert, del};
enum nucleotide {_A, _C, _G, _T, _X};

typedef vector< vector<double> > dmatrix;

inline bool isDelete(char x) {
  return x == DEL_CHAR;
}

inline bool isInsert(char x) {
     return islower(x);
}

inline bool isMatch(char x) {
  return isupper(x);
}

inline state_type getStateType(char x) {
  if (isMatch(x))
    return match;
  if (isInsert(x))
    return insert;
  if (isDelete(x))
    return del;
  assert(false);
  return match;    // To make the comiler stop whining.  
}

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
  }
  return _X;
}

bool isEmpty(const string& s) {
  for (string::const_iterator i=s.begin(); i!=s.end(); i++) {
    char c = *i;
    if (!(c == ' ' || c == '\t' || c == '\r' || c == '\n'))
      return false;
  }
  return true;
}

// Read data in data.  Each element of the vector corresponds to one line, and is broken into tokens by spaces.
// Skip blank lines and lines beginning with #.
void read_data(string file_name, vector<vector<string > >& lines){
  lines.clear();
  string temp_line;
  ifstream infile;
  infile.open (file_name.c_str());
  while(!infile.eof()) { // To get all the lines to the vector    
    getline(infile,temp_line); // Saves the line in temp_line.
    if (temp_line[0] != '#' && !isEmpty(temp_line)) {
      vector<string> one_line;   // to save the substrings in one single line
      istringstream iss(temp_line);
      while(iss){                // Split the line and push the substr to the vector
	string sub;
	iss >> sub;
	one_line.push_back(sub);
      }
      lines.push_back(one_line);
    }
  }
  infile.close();
}

void get_trans_count(vector< vector<string> >& lines, dmatrix& trans_count){
  trans_count = vector<vector<double> >(3, vector<double>(3, 0));
  vector< vector<string> >::iterator iter;

  for(iter=lines.begin()+1;iter!=lines.end();iter++){     // The "+1" is avoids processing the first line (which contains the ancestor)
    string temp_seq = (*iter)[ancestor_start];
    int i = 0;

    // KARRO: Modified the following for purposes of readability.  (No corrections to actual algorithm.)
    for (string::const_iterator s=temp_seq.begin(); s!=temp_seq.end()-1; i++) {
      state_type current_state = getStateType(*s);
      state_type next_state = getStateType(*(s+1));
      trans_count[current_state][next_state]++;
    }
  }

// KARRO: I'm removing the next section.  We will return the transition counts here, and calculate the probaibilities elsewhere.
//   dmatrix temp_trans_prob = trans_prob;    // to store the original matrix and help calculate the probs
//   for (int j=0; j<3; j++){
//     for(int k=0; k<3; k++){
//       trans_prob[j][k] = trans_prob[j][k]/(temp_trans_prob[j][0]+temp_trans_prob[j][1]+temp_trans_prob[j][2]);
// 	  //cout<<trans_prob[j][k]<<"\t";
//     }
//     //cout<<endl;
//   }
}

void better_emit_probs(vector< vector<string> > lines, vector<double> &ancestor_count, dmatrix &state_count, dmatrix &match_count, dmatrix &insert_count){   // both match and insert emit_matrix
  int num_columns = lines[0][1].length();             // the length of the ancestor
  string ancestor = lines[0][1];
 
  for(int i = 0; i < num_columns; i++){               // use the ancestor sequence to calculate P(x) vector
    nucleotide b = getNucleotide(ancestor[i]);
    if (b != _X) // Skipping degenerate bases
      ancestor_count[b]++;;
  }
  
  vector<int> num_each_col(num_columns, 0);           // to store the number of letters in each column of the ancestor, including "-" and lower cases
  for(vector<vector<string> >::iterator iter=lines.begin()+1;iter!=lines.end();iter++){
    int start_posi = atoi((*iter)[ancestor_start].c_str());
    string temp_seq = (*iter)[modern_sequence];
    int i = 0;
    int posi_inseq = 0;                            // the current position in this sequence
    
    while(i < temp_seq.length()){
      char b = temp_seq[i];
      nucleotide n = getNucleotide(b);
      if (isMatch(b) && n != _X) {
	state_count[start_posi + posi_inseq][0]++;
	match_count[start_posi + posi_inseq][n]++;
	posi_inseq++;
      }
      else if (isMatch(b) && n != _X) {
	state_count[start_posi + posi_inseq][1]++;
	insert_count[start_posi + posi_inseq][n]++;
      }
      else if (isDelete(b)) {
	state_count[start_posi + posi_inseq][2]++;
	posi_inseq++;
      }
      i++;
    }
  }
}

int main(int argc, char* argv[]) {

  string file;
  if (argc == 1) 
    file = "test.mal";
  else
    file = argv[1];

  // read_data
  vector< vector<string> > whole_file;
  read_data(file, whole_file);                    // read the file and store it into whole_file as strings

  string ancestor = whole_file[0][1];                  // to store the ancestor sequence
  int num_columns = ancestor.length();



  // Count samples
  dmatrix trans_count;   // Count number of transition types

  vector<double> ancestor_count(4, 0);   // Base counting matrix for the ancestor
  dmatrix state_count(num_columns, vector<double>(3, 0));           // to store the match/insert/delete matrix
  dmatrix match_count(num_columns, vector<double>(4, 0));          // to store the P(M/X) matrix
  dmatrix insert_count(num_columns, vector<double>(4, 0));

  dmatrix m_emit_prob(num_columns, vector<double>(4, 0));   // to store the match emit_prob matrix
  dmatrix i_emit_prob(num_columns, vector<double>(4, 0));   // to store the insert emit_prob matrix


  get_trans_count(whole_file, trans_count);    // transition prob matrix
  better_emit_probs(whole_file, ancestor_count, state_count, match_count, insert_count);      // pass by reference


  // Print Results
  cout<<ancestor<<endl;

  for(int i = 0; i < num_columns; i++){
	  for(int j = 0; j < 4; j++){
		  m_emit_prob[i][j] = match_count[i][j]/state_count[i][0];
		  i_emit_prob[i][j] = insert_count[i][j]/state_count[i][1];
		  cout<< m_emit_prob[i][j] <<"\t";
	  }
	  cout<<endl;
  }
  cout<<endl;

  /*pair <dmatrix, dmatrix> Pair_emit_matrices = better_emit_probs(whole_file, P_x, P_M, P_MX, P_IX);    
  dmatrix m_emit_prob = Pair_emit_matrices.first;      // match emission prob matrix
  dmatrix i_emit_prob = Pair_emit_matrices.second;     // insert emission prob matrix
  

  pHMM case1;
  case1.setName(file);
  case1.setAncestor(ancestor);
  for(int i = 0; i < 3; i++){          // set trans_probs
	  for(int j = 0; j < 3; j++){
		  case1.setTransProb(i, j, trans_prob[i][j]);
	  }
  }
  for(int i = 0; i < (int)m_emit_prob.size(); i++){
	  for(int j = 0; j < 4; j++){
		  case1.setEmissionProb(i, j, mat, m_emit_prob[i][j]);
		  case1.setEmissionProb(i, j, ins, i_emit_prob[i][j]);
	  }
  }*/
  
}

