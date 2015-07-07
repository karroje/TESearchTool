// input: ancestor seq; hashtable of chromosome; chromosome
// output: the chunks in the chromosome that has length of (anchor length + twice ancestor length)

#include "hunter.hh"
#include "sequence.hh"
#include <string>
#include <list>
#include <utility>
#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

vector<vector<string> > gener_chunks(string hash_file, string ances_file, string chr_file, string ances_name){
  //const char *ch_in = chr_file.c_str();
  //SEQ* s = seq_open(ch_in);    // store the chromosome into a string chr_str
  //seq_read(s);
    Sequence s(chr_file);
	Sequence ances_s(ances_file);
	while(ances_s.header().find(ances_name) == string::npos){
		ances_s.read();
	}

	//string chromosome_str = p;
	
	Hunter H(hash_file);            // hashtalble object
	int anchor_len = H.getWidth();
	vector<vector<string> > chunks;          // store all the output chunks here
	list<pair<int, int> > chunk_coords;         // store the coordinates (start, end) in the chromosome of each chunk, and use it to combine overlap

	for(int i=0; i<=((int)ances_s.size()-anchor_len); i++){       // for each anchor length seq in the ancestor
		// string temp_anchor = ances_seq.substr(i, anchor_len);  // store the anchor
		string temp_anchor = ances_s.subsequence(i, anchor_len);
		list<int> positions = H.getPositions(temp_anchor);     // store the positions of the anchor in the chromosome
		assert(H.getWidth() == (int)temp_anchor.size());
		for (list<int>::iterator p = positions.begin(); p != positions.end(); p++){   // for each position, generate the chunk
			int start_coor = 0;
			int chunk_len = anchor_len + (int)ances_s.size()*2;
			int end_coor = chunk_len;
			if(*p > chunk_len/2)
				start_coor = *p - chunk_len/2;
				end_coor = start_coor + chunk_len;
			//to reduce overlap by combining the overlapped coordinates
			bool overlap = false;
			for(list<pair<int, int> >::iterator p2 = chunk_coords.begin(); p2 != chunk_coords.end(); p2++){
				if((start_coor>=(*p2).first&&start_coor<=(*p2).second)||((*p2).first>=start_coor&&(*p2).first<=end_coor)){  //if chunks overlap
					(*p2).first = min((*p2).first, min((*p2).second, min(start_coor, end_coor)));
					(*p2).second = max((*p2).first, max((*p2).second, max(start_coor, end_coor)));
					overlap = true;
					break;
				}
			}
			if(!overlap){
				pair<int, int> new_chunk(start_coor, end_coor); 
				chunk_coords.push_back(new_chunk);
			}
		}
	}
	// To generate the chunks based on the coordinates
	for(list<pair<int, int> >::iterator p3 = chunk_coords.begin(); p3 != chunk_coords.end(); p3++){
		string temp_chunk = s.subsequence((*p3).first, anchor_len+(int)ances_s.size()*2);
		vector<string> one_chunk;
		one_chunk.push_back(temp_chunk);
		one_chunk.push_back(convertInt((*p3).first));
		one_chunk.push_back(convertInt(((*p3).first)+anchor_len+(int)ances_s.size()*2));
		one_chunk.push_back(convertInt((((*p3).first)+anchor_len+(int)ances_s.size()*2)-(*p3).first));
		chunks.push_back(one_chunk);
	}
	return chunks;
}

void gener_output(vector<vector<string> > chunks, string proc_id){
	ofstream myfile;
	string filename = proc_id + "_chunks.fa";
	myfile.open (filename.c_str());
	for (int i=0; i<chunks.size(); i++){
		//myfile << ">chunk" << "\t" << chunks[i][1] << "\t" << chunks[i][2] << "\n"; 
		myfile << ">chunk" << i << "\t" << chunks[i][1] << "\t" << chunks[i][2] << "\t" << chunks[i][3]<< "\n"; 
		myfile << chunks[i][0] << "\n";
	}
	myfile.close();

}

int main(int argc, char** argv) {
	string hash_file = argv[1];
	string ances_file = argv[2];
	string chr_file = argv[3];      // chromosome file
	string ances_name = argv[4];
	string proc_id = argv[5];

	vector<vector<string> > chunks = gener_chunks(hash_file, ances_file, chr_file, ances_name);
	gener_output(chunks, proc_id);

	//list<int> L = H.getPositions(ances_seq);
	//assert(H.getWidth() == ances_seq.size());
	
	//for (list<int>::iterator p = L.begin(); p != L.end(); p++)
	//	cout << *p << endl;
		
	return 0;
}
