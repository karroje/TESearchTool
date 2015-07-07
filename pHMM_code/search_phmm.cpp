/*
  @author Sean Wilkerson <wilkersm@muohio.edu>

*/


#include <iostream>
#include <fstream>
#include "matrix.hh"
#include "pHMM_b.hh"
#include <vector>
#include <string>
#include <cstring>
#include <cmath>
#include <map>
#include <limits>
//#include <boost/thread.hpp>

#include "hsp_hash/library.h"


using namespace std;

void* threadHandler(void* t);
vector<char> buildSequence(ifstream& i);

/*
   to make this where the multithreading happens,each thread needs its own pHMM instanstiation, as attempts to make
   the pHMM class threadsafe were a collassal failure, and while its possible, it will essentially throw away any 
   gained performance, 

   However if its infact feasible to make a new pHMM class for each thread we want to spawn then this is a nice place to seperate this
   for the final tool 

   
*/
/*
   For now arg 1 = filename to build to
   arg2 is filename to build from
   so usage is 

  search_hmm <filetosearch> <repeat_to_search_for> <anchor_start_coord> <anchor_end_coord>



 */

int main(int argc, char *argv[]){
  // for testing puposes only
  if (argc <= 1){
    printf("usage is  search_hmm <filetosearch> <repeat_to_search_for> <anchor_start_coord> <anchor_end_coord> <ancestor start>>\n");    exit(-1);
  }
  // TODO insert some getopt stuff in here
  
  //temp till above is done

  char * r_name = argv[2];
  char * gseq_in = argv[1];
  int astart = 0;//atoi(argv[5]); // ancestor linked coords
  int aend = 0;//atoi(argv[6]);// ancestor linked coords
  int anchstart = atoi(argv[3]);// anchor coords
  int anchend = atoi(argv[4]);// anchor coords
  
  
  
  cout << astart << endl;
  char h_name[256];
  char info_name[256];
  strcpy (h_name,r_name);
  strcpy (info_name,r_name);

  strcat(h_name, ".hmm");
  strcat(info_name, ".hmm.info");

  cout << "sequnce file = " << gseq_in << endl;
  cout << "hmm file = " << h_name << endl;
  cout << "hmm info = " << info_name << endl;

  ifstream infile (gseq_in, ifstream::in);
  ifstream hmmfile (h_name, ifstream::in);
  ifstream infofile (info_name, ifstream::in);

  //build ancestor table from

  if(! ( infile.good() && infile.is_open() && hmmfile.good() && hmmfile.is_open() && 
	 infofile.good() && infofile.is_open() ) ){
    cerr << "corrupt/non existant input file" << endl;
      return 0;
    }
  // note the program currently trusts the data to be valid and in the correct format
  // and only handles errors properly for missing/unopenable files
  // dump input file data to mem
  vector<char> observed = buildSequence(infile);

  
  pHMM model(hmmfile,infofile);
  
  /* The below sections needs reworked to work with model changes */

  
  // vector<char> p;

  // //search_pair forward (model.bDforward(observed,anchstart,anchend,astart,aend));
  // search_pair viterbi =  (model.bDviterbi(p,observed,anchstart,anchend,astart,aend));

  // // cout << "Foward Results (Log Space) :  " << forward.prob << " @" << forward.start << "-" << forward.end << endl; 
  // cout << "Viterbi Results (Log Space) :  " << viterbi.prob  << " @ " << viterbi.start << "-" << viterbi.end << endl; 
  // cout << "Viterbi Path: " << endl;
  // for (unsigned int i = 0; i < p.size(); ++i){
  //   cout << p[i];
  // }
  // cout << endl;
}


void* threadHandler(void *t){
  //currently unused, threading to be added after HSP database integration
  return NULL;
}


vector<char> buildSequence(ifstream& i){
  string linein;
  getline(i,linein);
  cout << linein << endl;
  
  int len = linein.length();
  vector<char> ret;

  for (int c = 0; c < len; ++c){
    ret.push_back(linein[c]);
  }
  
  return ret;
}
