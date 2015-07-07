#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctype.h>
#include <vector>
#include <list>
#include <math.h>

using namespace std;
#include "hunter.hh"
#include "seq.hh"
#define MAX_LINE_LENGTH 1048576

int main(int argc, char *argv[])
{    
    if(argc!=4)
    {
         cerr<<"Please give input file, output file and length of sequence!"<<endl;
         exit(1);
    }

    string input = argv[1];
    string output = argv[2];
    int len = atoi(argv[3]);
    
    const char *ch_in = input.c_str();
    const char *ch_out = output.c_str();
    
    int vec_len = 1 << len*2;
    int filter = vec_len - 1;
    vector<list<int> > vt(vec_len);

    SEQ* s = seq_open(ch_in); 
    seq_read(s);
    uchar* p = s->seq;
    
    int current_position = -1; //record the current position of the letter
    int code = 0; //unique code for each sequence
    int number; // unique code for one letter in the sequence
    int read_char = 0; //the number of continous non - N letters 
  
    while((char)*p!='\0')
    {
      current_position++;
                  
      read_char++;
                  
      char ch = (char)*p;
                  
      number = translation(ch);
                  
      if (number == 4)
	{
	  code = 0;
	  read_char = 0; 
	  p++;
	  continue;
	}
                  
      else 
	{
	  code = ((code<<2)& filter) | number;
	  p++;
	  if(read_char>=len)
	    {
	      vt[code].push_back(current_position-len+1);
	    }
	}
    }
  
    seq_close(s);
    ofstream fout(ch_out,ios::out | ios::binary);
  
    int list_size;
    //int char_n;
    int start = 0; 
    int int_size = sizeof(int);
    int element; // element in a list
    
    int reserve = int_size * 2; // reserved for sequence length and model
    int model = vec_len - 1; 
    // the current model is 11...11 which has the same length of the input sequence length
    
    int ele_start = int_size*2*vec_len + reserve; // spaces for the reserved spaces and the header 
    
    fout.write((char*)&len, int_size);
    fout.write((char*)&model, int_size);
    
/* output the result without code */
    for (int n = 0; n<vec_len; n++)
    {
        fout.seekp(int_size*2*n+reserve,ios::beg);
        list_size = vt[n].size();
        
        fout.write((char*)&list_size, int_size);
        fout.write((char*)&start, int_size);   
        fout.seekp(ele_start+start*int_size,ios::beg);
                 
        list<int>::iterator i = vt[n].begin();
        
        for(i = vt[n].begin(); i !=vt[n].end(); i++)
        {
            element=*i;
            fout.write((char*)&element, int_size);
        }
            
        start = start + list_size; 
    }
  

    fout.flush();
    fout.close();

    return 0;
}

