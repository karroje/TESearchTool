// Sample run command: -p test2.d/test.hmm -i test2.d/test.hmm.info -q test2.d/test.fa -h test2.d/test.hash -A 0.4 -C 0.1 -G 0.1 -T 0.4
#include <iostream>
#include <fstream>
#include "matrix.hh"
#include "pHMM.hh"
#include <vector>
#include <string>
#include <map>

#ifndef BUILD
#define BUILD

using namespace std;
// test

/*
  @author Sean Wilkerson <wilkersm@muohio.edu>

*/

/*
   For now arg 1 = filename to build to
   arg2 is filename to build from
   so usage is 

   build_hmm <destination> <source>
 */



int main(int argc, char *argv[]){
  // TODO insert some getopt stuff in here
  
  if (argc <= 1){
    printf("proper use is build_hmm <destination> <source>\n");
    exit(-1);
  }
  //temp till above is done

  char* fname_in = argv[2];
  char* fname_out = argv[1];
  
 
  ifstream infile (fname_in, ifstream::in);


  //build ancestor table from

  if(! ( infile.good() && infile.is_open() ) ){
    cerr << "corrupt/non existant input file" << endl;
      return 0;
  }
  // note the program currently trusts the data to be valid and in the correct format
  // and only handles errors properly for missing/unopenable files
  // dump input file data to mem

  string s;

  vector<string> data;
  while (!infile.eof() ){
    getline(infile,s);
    data.push_back(s);
  }

  pHMM newModel;
  string modelname(fname_out);
  newModel.buildFromAncestor(data,modelname);
  cout << "Done creating preliminary model, outputting to " << fname_out;
  cout << ".hmm and " << fname_out << ".hmm.info" << endl;

  if (newModel.printToFile(modelname)){
    cout << "Output Successful" << endl;
  }
  else
    cerr << "Error creating the output files" << endl;
  

}

#endif
