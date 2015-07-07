#include "pHMM.hh"
#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include "constants.hh"
#include "dna.h"
#include "util_fun.h"

using namespace std;

#define DEL_CHAR '-'

enum column_offsets {chr_location, chr_start, chr_stop, ancestor_start, ancestor_stop, modern_sequence};

typedef vector< vector<double> > dmatrix;
typedef vector<vector<int> > imatrix;

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
    return mat;
  if (isInsert(x))
    return ins;
  if (isDelete(x))
    return del;
  assert(false);
  return mat;    // To make the comiler stop whining.  
}

bool isEmpty(const string& s) {
  for (string::const_iterator i=s.begin(); i!=s.end(); i++) {
    char c = *i;
    if (!(c == ' ' || c == '\t' || c == '\r' || c == '\n'))
      return false;
  }
  return true;
}

void parseFileName(const string& input_name, string& file_dir, string& file_base) {
  int beginning = input_name.find_last_of('/');
  beginning = (beginning == (int)string::npos) ? 0 : beginning+1;

  int end = input_name.find_last_of('.');

  if (end == (int)string::npos || input_name.substr(end, input_name.size() - end) != ".mal") {
    cerr << "Input file does not appear to be a .mal file (based on name extension)." << endl;
    exit(1);
  }
  
  file_dir = input_name.substr(0, beginning);
  file_base = input_name.substr(beginning, end-beginning);
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
      vector<string> one_line;   // to save the substrs in one single line
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

void get_trans_count(vector< vector<string> >& lines, imatrix& trans_count){
  trans_count = vector<vector<int> >(3, vector<int>(3, 0));
  vector< vector<string> >::iterator iter;

  for(iter=lines.begin()+1;iter!=lines.end();iter++){     // The "+1" is avoids processing the first line (which contains the ancestor)
    string temp_seq = (*iter)[modern_sequence];

    // KARRO: Modified the following for purposes of readability.  (No corrections to actual algorithm.)
    for (string::const_iterator s=temp_seq.begin(); s!=temp_seq.end()-1; s++) {
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

void better_emit_probs(vector< vector<string> > lines, vector<int> &ancestor_count, imatrix &state_count, imatrix &match_count, imatrix &insert_count){   // both match and insert emit_matrix
  int num_columns = lines[0][1].length();             // the length of the ancestor
  string ancestor = lines[0][1];
 
  for(int i = 0; i < num_columns; i++){               // use the ancestor sequence to calculate P(x) vector
    nucleotide b = getNucleotide(ancestor[i]);
    if (b <= _T) // Skipping degenerate bases
      ancestor_count[b]++;;
  }
  
  vector<int> num_each_col(num_columns, 0);           // to store the number of letters in each column of the ancestor, including "-" and lower cases
  for(vector<vector<string> >::iterator iter=lines.begin()+1;iter!=lines.end();iter++){
    int start_posi = atoi((*iter)[ancestor_start].c_str());
    string temp_seq = (*iter)[modern_sequence];
    unsigned i = 0;
    int posi_inseq = -1;                            // the current position in this sequence
    
    while(i < temp_seq.length()){
      char b = temp_seq[i];
      nucleotide n = getNucleotide(b);
      if (isMatch(b) && n <= _T) {
	posi_inseq++;
	state_count[start_posi + posi_inseq][0]++;
	match_count[start_posi + posi_inseq][n]++;
      }
      else if (isInsert(b) && n <= _T) {
	state_count[start_posi + posi_inseq][1]++;
	insert_count[start_posi + posi_inseq][n]++;
      }
      else if (isDelete(b)) {
	posi_inseq++;
	state_count[start_posi + posi_inseq][2]++;
      }
      i++;
    }
  }
}

string help_string = "\
build_phmm:\n\
* Optional paramters:\n\
  -n <string>: Set name of hmm (default = input file name - .mal)\n\
  -o <string>: Output file (default = input file name - .mal + .hmm)\n\
* Required parameters:\n\
  1) Input file name\n\
";

int main(int argc, char* argv[]) {
  string hmm_name = "";
  string output_name = "";
  string file;

  int marker = 1;
  while (marker < argc && argv[marker][0] == '-') {
    string s = argv[marker];
    if (s == "--help") {
      cerr << help_string;
      exit(1);
    }
    if (s == "-n" || s == "-name")
      hmm_name = argv[++marker];
    if (s == "-o" || s == "-output")
      output_name = argv[++marker];
    ++marker;
  }

  if (marker == argc-1)
    file = argv[marker];
  else
    file = "test.mal";

  string file_dir;
  string file_base;
  parseFileName(file, file_dir, file_base);

  if (hmm_name == "")
    hmm_name = file_base;
  if (output_name == "")
    output_name = file_dir + file_base + ".hmm";

  // read_data
  vector< vector<string> > whole_file;
  read_data(file, whole_file);                    // read the file and store it into whole_file as strings

  string ancestor = whole_file[0][1];                  // to store the ancestor sequence
  int num_columns = ancestor.length();



  // Count samples
  imatrix trans_count;   // Count number of transition types

  vector<int> ancestor_count(4, 0);   // Base counting matrix for the ancestor
  imatrix state_count(num_columns, vector<int>(3, 0));           // to store the match/insert/delete matrix
  imatrix match_count(num_columns, vector<int>(4, 0));          // to store the P(M/X) matrix
  imatrix insert_count(num_columns, vector<int>(4, 0));

  imatrix m_emit_prob(num_columns, vector<int>(4, 0));   // to store the match emit_prob matrix
  imatrix i_emit_prob(num_columns, vector<int>(4, 0));   // to store the insert emit_prob matrix


  get_trans_count(whole_file, trans_count);    // transition prob matrix
  better_emit_probs(whole_file, ancestor_count, state_count, match_count, insert_count);      // pass by reference

  /***********
   * Create HMM
   ***********/
  pHMM H(hmm_name, ancestor.size());
  H.setAncestor(ancestor);

  // Set transition probabilities 
  for (int i=0; i < 3; i++) {
    double total = trans_count[i][0] + trans_count[i][1] + trans_count[i][2];
    for (int j=0; j < 3; j++)
      H.setTransProb(i, j, (total==0 ? neg_inf : log(trans_count[i][j] / total)));
  }

  // Set match emmission probabilities 
  for (int column=0; column < (int)ancestor.size(); column++) {
    double total = match_count[column][_A] + match_count[column][_C] + match_count[column][_G] + match_count[column][_T];
    for (int b = _A; b <= _T; b++)
      H.setEmissionProb(column, b, mat, total == 0 ? neg_inf : log(match_count[column][b] / total));
    H.setEmissionProb(column, _R, mat, logsum(H.emissionProb(column, _A, mat), H.emissionProb(column, _G, mat)));
    H.setEmissionProb(column, _Y, mat, logsum(H.emissionProb(column, _C, mat), H.emissionProb(column, _T, mat)));
    H.setEmissionProb(column, _N, mat, 0);
  }

  for (int column=0; column < (int)ancestor.size(); column++) {
    double total = insert_count[column][_A] + insert_count[column][_C] + insert_count[column][_G] + insert_count[column][_T];
    for (int b = _A; b <= _T; b++)
      H.setEmissionProb(column, b, ins, total == 0 ? neg_inf : log(insert_count[column][b] / total));
    H.setEmissionProb(column, _R, ins, logsum(H.emissionProb(column, _A, ins), H.emissionProb(column, _G, ins)));
    H.setEmissionProb(column, _Y, ins, logsum(H.emissionProb(column, _C, ins), H.emissionProb(column, _T, ins)));
    H.setEmissionProb(column, _N, ins, 0);
  }

  H.printToFile(output_name);
}

