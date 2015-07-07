
// input: ancestor seq; hashtable of chromosome; chromosome
// output: the chunks in the chromosome that has length of (anchor length + twice ancestor length)

//#include "hunter.hh"
//#include "sequence.hh"
#include <string>
#include <list>
#include <utility>
#include <vector>
#include <assert.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

int main(int argc, char** argv){
	vector <pair<int, int> > c1;
   vector <pair<int, int> >::iterator c1_Iter, temp_Iter;
   
   pair<int, int> temp_pair(10, 22);
   c1.push_back(temp_pair);
   pair<int, int> temp1_pair(20, 32);
   c1.push_back(temp1_pair);
   pair<int, int> temp2_pair(45, 57);
   c1.push_back(temp2_pair);
   pair<int, int> temp3_pair(46, 58);
   c1.push_back(temp3_pair);

   cout << "before erasing " << endl;
   for ( c1_Iter = c1.begin( ); c1_Iter != c1.end( ); c1_Iter++ )
	   cout << " " << c1_Iter->first << "   " << c1_Iter->second << endl;
   cout << endl;

   cout << "after conbining:" << endl;
   temp_Iter = c1.end() - 1;
   cout << temp_Iter->second <<endl;
   c1_Iter = c1.begin( );
   while ( c1_Iter->second != temp_Iter->second){
	   vector <pair<int, int> >::iterator c2_Iter;
	   c2_Iter = c1_Iter + 1;
	   if(c2_Iter->first - c1_Iter->second <= 12){     // if close enough, combine them
		   c1_Iter->second = c2_Iter->second;
		   c2_Iter = c1.erase(c2_Iter);
	   }
	   else{
		   c1_Iter++;
	   }
	   cout << "c1:  " <<c1_Iter->first << "   " << c1_Iter->second << "   c2:   "<<c2_Iter->first << "   " << c2_Iter->second << endl;
   }
   
   

   
}

/*
string convertInt(int number)
{
   stringstream ss;//create a stringstream
   ss << number;//add number to the stream
   return ss.str();//return a string with the contents of the stream
}

vector<vector<string> > gener_chunks(string hash_file, string ances_file, string chr_file, string ances_name, int len){
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
	vector<pair<int, int> > chunk_coords;         // store the coordinates (start, end) in the chromosome of each chunk, and use it to combine overlap
	vector<int> anchor_positions;            // store all the anchor positions
	for(int i=0; i<=((int)ances_s.size()-anchor_len); i++){       // for each anchor length seq in the ancestor
		// string temp_anchor = ances_seq.substr(i, anchor_len);  // store the anchor
		string temp_anchor = ances_s.subsequence(i, anchor_len);
		list<int> tmp_positions = H.getPositions(temp_anchor);     // store the positions of the anchor in the chromosome
		vector<int> positions(tmp_positions.begin(), tmp_positions.end());
		for (int j=0; j<positions.size(); j++)                     // push the coords of current anchor to the whole coords vector
			anchor_positions.push_back(positions[j]);
		//sort(positions.begin(), positions.end());
		assert(H.getWidth() == (int)temp_anchor.size());
	}
	sort(anchor_positions.begin(), anchor_positions.end());
	vector<pair<int, int> > anchor_coords;                         // store both the start and end coords of all the anchors
	for(int i=0; i<anchor_positions.size(); i++){
		pair<int, int> one_anchor(anchor_positions[i], anchor_positions[i]+anchor_len);
		anchor_coords.push_back(one_anchor);
	}



		//cout << "----------------------" << endl;
		for (vector<int>::iterator p = positions.begin(); p != positions.end(); p++){   // for each position, generate the chunk
			//cout << *p <<endl;       // for testing
			int start_coor = 0;
			int chunk_len = anchor_len + (len > 0 ? (int)ances_s.size()*len : -len);   // DEFINES CHUNK LENGTH
			int end_coor = chunk_len;
			if(*p > chunk_len/2){
				start_coor = *p - chunk_len/2;
				if(start_coor + chunk_len < (int)s.size())       
					end_coor = start_coor + chunk_len;
				else                                                    // if the chunk extends past the rear end of the genome, let the end of the genome be the end of the chunk
					end_coor = s.size() - 1;
			}
			// else: start_coor == 0, begining of chromosome


			// This section will reduce overlaping chunks by combining the overlapped coordinates into one chunk
			bool overlap = false;
			for(vector<pair<int, int> >::iterator p2 = chunk_coords.begin(); p2 != chunk_coords.end(); p2++){
				if((start_coor>p2->first && start_coor< p2->second) || (p2->first>start_coor && p2->first<end_coor)) {  //if chunks overlap
					p2->first = min(p2->first, min(p2->second, min(start_coor, end_coor)));
					p2->second = max(p2->first, max(p2->second, max(start_coor, end_coor)));
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
  for(vector<pair<int, int> >::iterator p3 = chunk_coords.begin(); p3 != chunk_coords.end(); p3++){
    cout << p3->first << '\t' << p3->second <<endl;
    string temp_chunk = s.subsequence(p3->first, anchor_len+(int)ances_s.size()*2);
    vector<string> one_chunk;
    one_chunk.push_back(temp_chunk);
    one_chunk.push_back(convertInt(p3->first));
    one_chunk.push_back(convertInt((p3->first)+anchor_len+(int)ances_s.size()*2));
    one_chunk.push_back(convertInt(((p3->first)+anchor_len+(int)ances_s.size()*2)-p3->first));
    chunks.push_back(one_chunk);
  }
  return chunks;
}

void gener_output(vector<vector<string> > chunks, string proc_id){
	ofstream myfile;
	string filename = proc_id + "_chunks.fa";
	myfile.open (filename.c_str());
	for (int i=0; i < (int)chunks.size(); i++){
		//myfile << ">chunk" << "\t" << chunks[i][1] << "\t" << chunks[i][2] << "\n"; 
		myfile << ">chunk" << i << "\t" << chunks[i][1] << "\t" << chunks[i][2] << "\t" << chunks[i][3]<< "\n"; 
		myfile << chunks[i][0] << "\n";
	}
	myfile.close();

}

int main(int argc, char** argv) {
  int len = 2;    // If positive, chunk length = len*ancestor_size.
                  // If negative, chunk length = ancestor_size + (-1)*len
  int marker = 1;
  while (marker < argc && argv[marker][0] == '-') {
    string s = argv[marker++];
    if (s == "-c" || s == "-cl" || s == "--chunk_len") {
      len = atoi(argv[marker++]);
    }
    else if (s == "-h" || s == "--hard_chunk_limit") {
      len = -1*atoi(argv[marker++]);
    }
  }

  if (marker + 5 != argc) {
    cerr << "chunks.cpp -- bad command line" << endl;
    exit(1);
  }
  string hash_file = argv[marker];
  string ances_file = argv[marker+1];
  string chr_file = argv[marker+2];      // chromosome file
  string ances_name = argv[marker+3];
  string proc_id = argv[marker+4];

  vector<vector<string> > chunks = gener_chunks(hash_file, ances_file, chr_file, ances_name, len);
  gener_output(chunks, proc_id);

  return 0;
}
*/









/*
void testTop(int t){
	--t;
}

int main(int argc, char** argv){

	list <int> c1;
   list <int>::iterator c1_Iter;
   vector <int>::iterator c2_Iter;
   
   c1.push_back( 20 );
   c1.push_back( 10 );
   c1.push_back( 30 );

   cout << "Before sorting: c1 =";
   for ( c1_Iter = c1.begin( ); c1_Iter != c1.end( ); c1_Iter++ )
      cout << " " << *c1_Iter;
   cout << endl;

   c1.sort( );
   cout << "After sorting c1 =";
   for ( c1_Iter = c1.begin( ); c1_Iter != c1.end( ); c1_Iter++ )
      cout << " " << *c1_Iter;
   cout << endl;

   vector <int> c2(c1.size());
   copy(c1.begin(), c1.end(), c2.begin());
   for ( c2_Iter = c2.begin( ); c2_Iter != c2.end( ); c2_Iter++ )
      cout << " " << *c2_Iter;
   cout << endl;
 */  
	/*int *a;
	a = new int[4];
	a[0] = 1;
	cout<<a[0]<<endl;

	hash_map<int, int> hm;
	
	int arr1[] = {1, 3, 4, 5};
	int arr2[] = {4, 6, 7, 8, 9};
	int arr3[] = {2, 4, 10};

	for(int i=0; i<4; i++)
		hm[arr1[i]] = 1;
	for(int i=0; i<5; i++){
		if(hm[arr2[i]] != 0)
			hm[arr2[i]]++;
		else
			hm[arr2[i]] = 1;
	}
	for(int i=0; i<3; i++){
		if(hm[arr3[i]] != 0)
			hm[arr3[i]]++;
		else
			hm[arr3[i]] = 1;
	}

	hash_map<int, int>::iterator iter;
	for(iter = hm.begin(); iter != hm.end(); iter++){
		cout<<"hmKey: "<<iter->first<<"   hmVal: "<<iter->second<<endl;
		if(iter->second == 3)
			cout<<"common value: "<<iter->first<<endl;
	}*/
	/*
	hm["wow"] = 5;
	hm["yeah"] = 8;
	cout<<"hm[\"1\"]: "<<hm["1"]<<endl;
	hash_map<string, int>::iterator iter;
	for(iter = hm.begin(); iter != hm.end(); iter++){
		cout<<"hmKey: "<<iter->first<<"   hmVal: "<<iter->second<<endl;
	}*/
	/*
	cout<<"Please enter: "<<endl;
	int a;
	cin>>a;

	for(int i=1; i<=a; i++){
		for(int j=1; j<=i; j++)
			cout<<j;
		cout<<endl;
	}
	cout<<"\n";
	cout<<"Please enter another:"<<endl;
	int b;
	cin>>b;
	for(int i=1; i<=b; i++){
		for(int j=1; j<=i; j++)
			cout<<"*";
		cout<<endl;
	}*/
	//int top = 5;
	//testTop(top);
	//cout<<top<<endl;
	/*
	string str;
	cin>>str;
	cout<<str+" _hello"<<endl;
	string filename = str+" _hello";
	ofstream myfile;
	myfile.open(filename);

	return 0;*/