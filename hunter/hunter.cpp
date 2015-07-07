#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <ctype.h>
#include <list>
#include <vector>
using namespace std;
#include "hunter.hh"

Triple::Triple(int f, int s, int c)
{
   first_s = f;         //Starting pos of first sequence  
   second_s = s;        //Starting pos of second sequence
   common_length = c;   //Length of the alignment.  (No insertions of deletions)
}


ostream& operator<<(ostream& out, const Triple& T) {
  out << "(" << T.first_s << "," << T.second_s << "," << T.common_length << ")";
  return out;
}

/*
fn is the name of the fasta file
w is the width of the anchors
 */
Hunter::Hunter(string fn){
     filename = fn;
     fp.open(filename.c_str(), ios::in | ios::binary);

     fp.read((char *)&width, int_size);
     fp.read((char *)&model, int_size);
     
     // vec_len = (int)pow(2,width*2);
     vec_len = 1 << width*2;
     reserve = int_size << 1;
     header_space = int_size2*vec_len + reserve;
     }
     
int Hunter::getWidth()const
{
     return width;
}

int Hunter::getModel()const
{
     return model;
}

list<int> Hunter::getPositions(string seq)
{
     int seq_length = seq.length();
     int code = 0;
     int list_len;//length of the list
     int start_position; // start position of the list in binary file
     int element;// element in the list
     list<int> position_list;
     
     if (seq_length!=width)
     {
         cerr<<"The length of the sequence is wrong!"<<endl;
         exit(1);
     }                 
      

     code = translate_seq(seq.c_str());

     fp.seekg(code*int_size2+reserve,ios::beg);

     fp.read((char *)&list_len, int_size);
     fp.read((char *)&start_position, int_size);

     fp.seekg(header_space+start_position*int_size);

     for (int i = 0; i<list_len; i++)
     {
         fp.read((char*)&element, int_size); 
         position_list.push_back(element);
     }
     
     
     /*
     list<int>::iterator i = position_list.begin();
        
     for(i = position_list.begin(); i !=position_list.end(); i++)
     {
        element=*i;
        cout<<element<<" ";
     }
     */
     return position_list;
}

//list<Triple> Hunter::getTriples(string s, int l){
list<Triple> Hunter::getTriples(string s){
  int str_len = s.length();
  vector<list<int> > vt(str_len-width+1);
  string sub_str;
  
  list<int> test_list;
    
  for(int i = 0; i < str_len-width+1;i++)
    {
      sub_str = s.substr(i,width);
      //std::cout<<endl<<"The substring is: "<< sub_str<<":";
      //vt[i]=getPositions(sub_str);
      test_list = getPositions(sub_str);
      if (test_list.empty())
          vt[i].push_back(-25);
      else vt[i] = test_list;      
    }
    
  int vector_size = vt.size();
  //std::cout<<endl<<"The size of the vector is: "<<vector_size<<endl;
  
  vector<Triple> position_vector[2];
  
  vector<list<Triple> > result_vector(str_len-width+1);
  
  int mark = 0;
  
  //std::cout<<endl<<"Initializing..."<<endl;
  for(list<int>::iterator m = vt[0].begin();m!=vt[0].end();m++)
  {
      int first = *m;
      //cout<<endl<<"put {"<<first<<",0,1} into vector";
      Triple t = Triple(first,0,1);
      position_vector[0].push_back(t);
  }
  
  //put array into the vector
  //std::cout<<endl<<"Start to put triples into vector..."<<endl;
  for(int i=1;i<vector_size;i++)
  {
     //std::cout<<"i="<<i<<endl;
     position_vector[1-mark].clear();
     
     list<int>::iterator m = vt[i].begin();
     
     vector<Triple>::iterator it = position_vector[mark].begin();
     
     while((m!=vt[i].end())||(it!=position_vector[mark].end()))
     {
       
       int first = *m;
       
       int common_len = 1;

       Triple tr = *it;
       
       if(((tr.first_s+1)==first)&&((tr.second_s+1)==i))// find a continued element
       {
           
           common_len=tr.common_length+1;
           Triple tri = Triple(first, i, common_len);
           //cout<<endl<<"push_back triple 1("<<tri.first_s<<","<<tri.second_s<<","<<tri.common_length<<") "<<endl;
           position_vector[1-mark].push_back(tri);
          // Triple result_tr = Triple(tr.first_s-tr.common_length+1,tr.second_s-tr.common_length+1,tr.common_length+width-1);
          // result_vector[tr.common_length-1].push_back(result_tr);
           m++;
           it++;
           
       }
       
       else if (((it!=position_vector[mark].end())&&(tr.first_s<=first))||(m==vt[i].end()))//
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
                   Triple result_tr = Triple(tr.first_s-tr.common_length+1,tr.second_s-tr.common_length+1,tr.common_length+width-1);
                   if (result_tr.first_s != -25)
                       result_vector[tr.common_length-1].push_back(result_tr);
                       //cout<<endl<<"push_back triple 2("<<result_tr.first_s<<","<<result_tr.second_s<<","<<result_tr.common_length<<") "<<endl;
                   it++;
               }
           }
           else
           {
               Triple result_tr = Triple(tr.first_s-tr.common_length+1,tr.second_s-tr.common_length+1,tr.common_length+width-1);
               if (result_tr.first_s != -25)
                   result_vector[tr.common_length-1].push_back(result_tr);
               //cout<<endl<<"push_back triple 2("<<result_tr.first_s<<","<<result_tr.second_s<<","<<result_tr.common_length<<") "<<endl;
               it++;
           }
       }
       
       else if (((tr.first_s>first)&&(m!=vt[i].end()))||(it==position_vector[mark].end()))
       {
            
            Triple tri = Triple(first, i, common_len);
            position_vector[1-mark].push_back(tri);
            //cout<<endl<<"push_back triple 3("<<tri.first_s<<","<<tri.second_s<<","<<tri.common_length<<") "<<endl;
            m++;
       }    
      
      else m++; 
      //cout<<"while end"<<endl;
     }
     
     mark = 1-mark;
     
     if(i==vector_size-1)
     {
         //std::cout<<"push_back last vector"<<endl;
         for(vector<Triple>::iterator it = position_vector[mark].begin();it!=position_vector[mark].end();it++)  
         {
            Triple trip=*it;
            Triple result_tr = Triple(trip.first_s-trip.common_length+1,trip.second_s-trip.common_length+1,trip.common_length+width-1);
            //std::cout<<"("<<result_tr.first_s<<","<<result_tr.second_s<<","<<result_tr.common_length<<") ";
            if (result_tr.first_s != -25)
                result_vector[trip.common_length-1].push_back(result_tr);
         }
     }             
  } 
  
/*  std::cout<<endl<<"Result Vector..."<<endl;
  int result_vector_size = result_vector.size();
  for(int i=0;i<result_vector_size;i++)
  {
      std::cout<<endl<<"Result Vector["<<i<<"] is ";
      for(list<Triple>::iterator ti = result_vector[i].begin();ti!=result_vector[i].end();ti++)
      {
          Triple result_triple = *ti;
          //if (result_triple.first_s == -25)
              //continue;
          //else 
              std::cout<<"("<<result_triple.first_s<<","<<result_triple.second_s<<","<<result_triple.common_length<<") ";
      }
   }
   
   return result_vector[l-width];
*/

  list<Triple> result_list;
  //std::cout<<endl<<"Result Vector..."<<endl;
  int result_vector_size = result_vector.size();
  for(int i=0;i<result_vector_size;i++)
  {
      //std::cout<<endl<<"Result Vector["<<i<<"] is ";
      for(list<Triple>::iterator ti = result_vector[i].begin();ti!=result_vector[i].end();ti++)
      {
          Triple result_triple = *ti;
          result_list.push_back(result_triple);
          //std::cout<<"("<<result_triple.first_s<<","<<result_triple.second_s<<","<<result_triple.common_length<<") ";
      }
   }
   
  // return result_vector[l-width];
  return result_list;
}              

/*int main()
{
    return 0;
}*/

