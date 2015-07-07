#include "hunter.hh"
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctype.h>
#include <list>
#include <vector>
#include <string.h>
#include "HashTable.h"
//#include <iostream>

using namespace std;

Triple::Triple(int f, int s, int c)
{
  first_s = f; //Starting pos of first sequence  
  second_s = s;//Starting pos of second sequence
  common_length = c;//Length of the alignment.  (No insertions of deletions)
}


ostream& operator<<(ostream& out, const Triple& T) {
  out << "(" << T.first_s << "," << T.second_s << "," << T.common_length << ")";
  return out;
}
/*
fn is the name of the fasta file
w is the width of the anchors
 */
Hunter::Hunter(string fn,int w){
     filename = fn;
     width = w;
     buildHash();
     handle = hashes.begin();
}

void Hunter::buildHash(){
  //Should I keep the fragments in memory?
  list<string> frags = readFasta(filename);
  HashTable* h = new HashTable();
  for(list<string>::const_iterator it = frags.begin();
      it != frags.end(); it++){
    for(int i = 0; i<(*it).length()-width;i++){
      h->insert((*it).substring(i,i+width),i);
    }
    hashes.push_back(h);
  }
  delete frags;
}

list<string> Hunter::readFasta(string seqfile){
  ifstream fin((char*)seqfile.c_str());
  string tmp;
  string cur_seq = "";
  list<string> seqs;
  getline(fin,tmp);
  while(!fin.eof()){
    while(tmp[0]!=">"){
      cur_seq+=tmp;
      getline(fin,tmp);
    }
    getline(fin,tmp);
    seqs.push_back(cur_seq);
    cur_seq = "";
  }
  return seqs;
}


int Hunter::getWidth()const
{
     return width;
}

int Hunter::getModel()const
{
     return model;
}


list<Triple> Hunter::getTriples(string s){
  list<Triple> hits;
  string subseq;
  list<int> pos //Positions on the fragment
  for(int i = 0; i<s.length()-width;i++){
    subseq = s.substr(i,i+width);
    pos = hash_it->find(subseq);
    for(list<int>::const_iterator it = pos.begin(); 
	it != pos.end(); it++){
      hits.push_back(Triple(i,(*it),width));
    }
  }
  return hits;
  
  /*
  int str_len = s.length();
  vector<list<int> > vt(str_len-width+1); 
  string sub_str;
  
  list<int> test_list;
    
  for(int i = 0; i < str_len-width+1;i++)
    {
      sub_str = s.substr(i,width);
      test_list = getPositions(sub_str);
      if (test_list.empty())
          vt[i].push_back(-25);
      else vt[i] = test_list;      
    }
    
  int vector_size = vt.size();
  
  vector<Triple> position_vector[2];
  
  vector<list<Triple> > result_vector(str_len-width+1);
  
  int mark = 0;
  
  for(list<int>::iterator m = vt[0].begin();m!=vt[0].end();m++)
  {
      int first = *m;
      Triple t = Triple(first,0,1);
      position_vector[0].push_back(t);
  }
  
  for(int i=1;i<vector_size;i++)
  {
     position_vector[1-mark].clear();
     
     list<int>::iterator m = vt[i].begin();
     
     vector<Triple>::iterator it = position_vector[mark].begin();
     int it_flag = 0;//indicates if it==end()
     int vt_flag = 0;//indicates if vt==end()
     //JAMIE: Changed || to && below
     //while((m!=vt[i].end())||(it!=position_vector[mark].end()))
     while(it_flag!=1 and m!=vt[i].end())
     { 
       if(m==vt[i].end()){vt_flag=1;}
       if(it==position_vector[mark].end()){it_flag=1;}
        
       int first = *m;
       
       int common_len = 1;

       Triple tr = *it;
       
       if(((tr.first_s+1)==first)&&
	  ((tr.second_s+1)==i))// find a continued element
       {
           
           common_len=tr.common_length+1;
           Triple tri = Triple(first, i, common_len);
           position_vector[1-mark].push_back(tri);
	   m++;
           it++;
           
       }
       
       else if (((it!=position_vector[mark].end())&&
		 (tr.first_s<=first))||
		(m==vt[i].end()))//
       {
           if (tr.first_s == first)
           {
               m++;
               int next =*m;
               m--;
               if (next == tr.first_s+1)
               {
                   Triple tri = Triple( first, i, common_len );
                   position_vector[1-mark].push_back(tri);
                   m++;
               }
               else 
               {
                   Triple result_tr = Triple(tr.first_s-tr.common_length+1,
					     tr.second_s-tr.common_length+1,
					     tr.common_length+width-1);
                   if (result_tr.first_s != -25)
                       result_vector[tr.common_length-1].push_back(result_tr);
                   it++;
               }
           }
           else
           {
               Triple result_tr = Triple(tr.first_s-tr.common_length+1,
					 tr.second_s-tr.common_length+1,
					 tr.common_length+width-1);
               if (result_tr.first_s != -25 && first != -25)
                   result_vector[tr.common_length-1].push_back(result_tr);
               it++;
           }
       }
       
       else if (((tr.first_s>first)&&
		 (m!=vt[i].end()))||
		(it==position_vector[mark].end()))
       {
            
            Triple tri = Triple(first, i, common_len);
            position_vector[1-mark].push_back(tri);
            m++;
       }    
      
      else m++; 
     }
     vt_flag = 0;
     it_flag = 0;
     
     mark = 1-mark;
     
     if(i==vector_size-1)
     {
         for(vector<Triple>::iterator it = position_vector[mark].begin();
	     it!=position_vector[mark].end();
	     it++)  
         {
            Triple trip=*it;
            Triple result_tr = Triple(trip.first_s-trip.common_length+1,
				      trip.second_s-trip.common_length+1,
				      trip.common_length+width-1);
            if (result_tr.first_s != -25)
                result_vector[trip.common_length-1].push_back(result_tr);
         }
     }
     //cout<<"Flags: "<<m!=vt[i].end()<<" "<<it!=position_vector[mark].end();

  } 
  
  

  list<Triple> result_list;
  int result_vector_size = result_vector.size();
  for(int i=0;i<result_vector_size;i++)
  {
      for(list<Triple>::iterator ti = result_vector[i].begin();
	  ti!=result_vector[i].end();
	  ti++)
      {
          Triple result_triple = *ti;
          result_list.push_back(result_triple);
      }
   }
   
  return result_list;*/
}              

