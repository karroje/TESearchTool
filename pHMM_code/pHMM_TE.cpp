// Sample run: -p test2.d/test.hmm -i test2.d/test.hmm.info -q test2.d/test.fa -h test2.d/test.hash -A 0.5 -C 0.1 -G 0.1 -T 0.4

#include <vector>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cmath>
#include <string>
#include <sstream>
#include <utility>
#include <cstring>
#include <unistd.h>
#include <stdio.h>
#include <ctype.h>
#include <getopt.h>
#include <boost/tuple/tuple_io.hpp>

#include "hunter/hunter.hh"
#include "seq.hh"
#include "pHMM.hh"
#include "matrix.hh"
#include "utils.h"
#include "endProbFunc.hh"

/*
  command line usage, and common error strings
 */
using namespace std;
using namespace boost;

static char CL_ARGS[] = "\
Command line parameters:  \n\
   -p <file>   An HMM for a specific family (required). \n\
   -q <file>   A sequence file for a chromosome (in .fa or .nib format) (required). \n\
   -h <file>   A .hash file for the same chromosome. (required)\n\
   -s <start coordinate>  genomic start coordinate (optional) \n\
   -f <end coordinate> genomic end coordinate (optional).\n\
   -A/C/G/T <double> Set the background distribution for each base (optional, default = (1,1,1,1))\n\
   -B <axbxcxd> Set the background distributions at one time\n\
   -e <endprob type>  Assumed endpoint distribution (default = fixed)\n\
   -H displays this help text\n";

// end static strings


//*******************
// prameters
string hashfile; // -h
string seqfile; // -q
string hmmfile; // -p
int sflag = 0;
string coord_S; // -s
int fflag = 0;
string coord_F; // -f

string EPType = "fixed";
endProbFuncFactory EPFact;

double background_frequency_total;
double background_frequencies[4] = {1,1,1,1};
double normalized_bg_freq[4]; // Elements set with -A, -C, -G and -T

SEQ* sequence;
Hunter* hunt;
pHMM* model;
string ancestor;
int aLen;
vector<string> anchor_set;
list<Triple> hunterHits;
long start_C, end_C;

// void parse_cmdline(int argc, char** argv,
// 		   int& hflag, char*& hashfile,
// 		   int& qflag, char*& seqfile,
// 		   int& pflag, char*& hmmfile,
// 		   int& sflag, char*& coord_S,
// 		   int& fflag, char*& coord_F,
// 		   double background_frequencies[],
// 		   double normalized_bg_freq[],
// 		   string& EPType,
// 		   int& Hflag) {
void parse_cmdline(int argc, char** argv) {
    opterr = 0;

    // defaults
    hashfile = (string)"test.hash";
    hmmfile = "test.hmm";
    seqfile = "test.fa";
    EPType = "uniform";
    background_frequencies[0] = 0.1;
    background_frequencies[1] = 0.4;
    background_frequencies[2] = 0.4;    
    background_frequencies[3] = 0.1;


    int c;
    while ((c = getopt (argc, argv, "A:C:G:T:p:q:h:i:s:e:f:H:B:")) != -1) {
      switch (c)
	{
	case 'p':
	  hmmfile = optarg;
	  break;
	case 'q':
	  seqfile = optarg;
	  break;
	case 'h':
	  hashfile = optarg;
	  break;
	case 's':
	  sflag = 1;
	  coord_S = optarg;
	  break;
	case 'f':
	  fflag = 1;
	  coord_F = optarg;
	  break;
	case 'e':
	  EPType = optarg;
	  break;
	case 'B': // Set background distribution
	  background_frequencies[0] = atof(strtok(optarg, "x"));
	  background_frequencies[1] = atof(strtok(NULL, "x"));
	  background_frequencies[2] = atof(strtok(NULL, "x"));
	  background_frequencies[3] = atof(strtok(NULL, "\0"));
	  break;
	case 'A':
	  background_frequencies[0] = log(atof(optarg));
	  break;
	case 'C':
	  background_frequencies[1] = log(atof(optarg));
	  break;
	case 'G':
	  background_frequencies[2] = log(atof(optarg));
	  break;
	case 'T':
	  background_frequencies[3] = log(atof(optarg));
	  break;
	case 'H':
	  printf("%s", CL_ARGS);
	  exit(-1);
	  break;
	case '?':
	  if (optopt == 'p' || optopt == 'q' || optopt == 'h' || optopt == 's' || optopt == 'f')
	    fprintf (stderr, "Option -%c requires an argument.\n", optopt);
	  else if (isprint (optopt))
	    fprintf (stderr, "Unknown option `-%c'.\n", optopt);
	  else
	    fprintf (stderr,
		     "Unknown option character `\\x%x'.\n",
		     optopt);
	  exit(1);
	default:
	  printf("%s", CL_ARGS);
	  abort ();
	}
    }

    if(!(fflag == sflag)){
      die("Start and finish parameteres must be provided in pairs, run with -H for help\n");
    }
}		      

void setup() {
  EPFact.setType(EPType);
  background_frequency_total = 0;
  for (int i=0; i < 4; i++)
    background_frequency_total += background_frequencies[i];
  for (int i=0; i < 4; i++)
    normalized_bg_freq[i] = log(background_frequencies[i] / background_frequency_total);

  sequence = seq_open(seqfile.c_str());
  
  if (!seq_read(sequence)){
    die("error processing reading sequence file\n");
  }

  hunt = new Hunter(hashfile);
  model = new pHMM((char*)(hmmfile.c_str()));
  ancestor = model->getAncestorSequence();
  aLen = ancestor.size();

  for (int i = 0 ; i < aLen - hunt->getWidth(); ++i){
    anchor_set.push_back(ancestor.substr(i,hunt->getWidth()));
  }

  hunterHits = hunt->getTriples(ancestor);
  start_C = coord_S != "" ? atol(coord_S.c_str()) : 0;
  end_C = coord_F != "" ? atol(coord_F.c_str()) : aLen-1;
}


int main(int argc, char **argv){


    //***********************************************************************
    // Setup

    //parse_cmdline(argc, argv, hflag, hashfile, qflag, seqfile, pflag, hmmfile, sflag, coord_S, fflag, coord_F, background_frequencies, normalized_bg_freq, EPType, Hflag);
    parse_cmdline(argc, argv);
    setup();


    list<Triple>::iterator it;

    Alignment Best;
    tuple<int,int,double> BestF = make_tuple(0,0,neg_inf);

    list<Triple>::iterator it_best;
    list<Triple>::iterator it_best_f;
    for (it = hunterHits.begin(); it != hunterHits.end(); it++){
      /* If flags are not set or flags are set and the triple is inside those bounds (if
       the flags produce unrealisitc bounds this isn't checked but it will just rule everything out) */
      if (!fflag || ((it->first_s >= start_C ) && (( it->first_s + it->common_length) < end_C))){
	endProbFunc* EPFunc = EPFact.createEndProbFunc(model->length(), it->second_s, it->second_s + it->common_length - 1);

	tuple<int,int,double> F = model->bDforward(sequence, it->first_s, it->first_s + it->common_length - 1, it->second_s, it->second_s + it->common_length - 1, normalized_bg_freq, EPFunc, false);
	if (F.get<2>() > BestF.get<2>()) {
	  it_best_f = it;
	  BestF = F;
	}

	Alignment A = model->bDviterbi(sequence, it->first_s, it->first_s + it->common_length - 1, it->second_s, it->second_s + it->common_length - 1, normalized_bg_freq, EPFunc, false);
	//cout << it->first_s << " " << it->second_s << " " << it->common_length << " " << A.score << endl;

	if (A.score > Best.score) {
	  it_best = it;
	  Best = A;
	}
	delete EPFunc;
      }
    }

    cout << "Forward:\n";
    cout << "Anchor: " << it_best->first_s << " " << it_best->second_s << " " << it_best->common_length << endl;
    cout << "Solution: " << BestF << endl << endl;

    cout << "Viterbi:\n";
    cout << "Anchor: " << it_best->first_s << " " << it_best->second_s << " " << it_best->common_length << endl;
    pair<int,int>  m_ends = Best.modern_ends();
    pair<int,int>  a_ends = Best.ancestor_ends();
    cout << "Modern: " << make_tuple(m_ends.first, m_ends.second) << endl;
    cout << "Ancestor: " << make_tuple(a_ends.first, a_ends.second) << endl;

    cout << "Score: " << exp(Best.score) << endl;
    if (Best.score > log(0)) 
      Best.print_alignment(cout, sequence, model->getAncestorSequence());
    cout << "\n";


    // output
    /*  this is still up in the air  right now the output is the printf above which just
    dumps everything above a certain arbitrary threshold*/

 return 0;
}
