"""
@package docstring
The purpose of this program is to parse the output of the hmmsearch program from the HMMER package.  
The HMMparser class encases all of the hits from the hmmsearch program 
The HMMRecord class holds all of the information for a particular hit
The Domain class holds all of the information about potential alignments
Once the program is parsed, the user has the option to mask the adapter sequences
Writes the result to _hmmparse.out
""" 
import re
import sys

class HMMparser(object):
    """
    Parses the output of hmmsearch
    """
    def __init__(self,file_input):
        self.file_in = file_input
        self.num_records = 0
        self.hmm_records = []
        self.parseHMMoutput()
        self.num_records = len(self.hmm_records)
    def __str__(self):
        string = ""
        for i in range(0,len(self.hmm_records)):
            #print i
            string+=str(self.hmm_records[i])+'\n'
        return string
    def __del__(self):
        del self.hmm_records
    """
    This fills in the data structure with alignment information
    """
    def parseHMMoutput(self):
        handle = open(self.file_in,'r')
        lines = handle.readlines()
        #print lines
        self.build_parser(lines)
        index = self.num_records
        i = 0       
        #debug = open('debug_args.txt','w')
        #debug_parse_output = open('debug_parse_output.txt','w')
        #while(index<len(lines)-15):
        while(i<self.num_records):
            while([s for s in lines[index].split(' ') if s!=''][0]!='>>'):
                if index>len(lines)-15:
                    return
                index+=1
            #if([s for s in lines[index].split(' ') if s!=''][1]!=self.hmm_records[i].description):
               #print ">> " + self.hmm_records[i].description
               #print lines[index]+'\n'
            
            if([s for s in lines[index+1].split(' ') if s!=''][0]=='[No'):#Domain is nonexistant
                #del self.hmm_records[i]
                i+=1
                index+=1
                continue
                
            while([s for s in lines[index].split(' ') if s!=''][0].isdigit()==False):
                if index>len(lines)-15:
                    return
                index+=1
            for d in range(0,self.hmm_records[i].num_domains):#Initializes all domains
                d_args = lines[index].split(' ')
                d_args = [s for s in d_args if s!='!' and s!='[.' and s!='[]' and s!='.]'\
                          and s!='..' and s!='?' and s!='']
                #print index
                #if(len(d_args)<8):
                #    print d_args
                #    print index
                    
                dom = Domain(score = float(d_args[1]),
                             bias = float(d_args[2]),
                             cEvalue = float(d_args[3]),
                             iEvalue = float(d_args[4]),
                             query_start = int(d_args[7]),
                             query = '',
                             query_end = int(d_args[8]),
                             sbjct_start = int(d_args[5]),
                             sbjct = '',
                             sbjct_end = int(d_args[6]))
                self.hmm_records[i].domains.append(dom)
                index+=1
            for d in range(0,self.hmm_records[i].num_domains):
                while([s for s in lines[index].split(' ') if s!=''][0]!='=='):
                    index+=1
                
                index+=1
                sbjct_args = lines[index].split(' ')
                sbjct_args = [s for s in sbjct_args if s!='']
                query_args = lines[index+2].split(' ')
                query_args = [s for s in query_args if s!='']
                self.hmm_records[i].domains[d].sbjct = sbjct_args[2].upper()           
                self.hmm_records[i].domains[d].query = query_args[2].upper()
            #print>>debug_parse_output,"Index: "+str(index)+" Record: "+str(i) 
            #print>>debug_parse_output,self.hmm_records[i]
            i+=1
        
        #Find any faulty records in the data structure and remove them
        n = 0
        while(n<len(self.hmm_records)):
            if(len(self.hmm_records[n].domains)==0):
               self.hmm_records.pop(n)
            n+=1
            #print n
        self.num_record = len(self.hmm_records)
        del lines
    """
    This builds the data structure that holds a list of HMMRecords
    """
    def build_parser(self,lines):
        #debug_build = open('debug_build_parser.txt','w')
        index = 0
        #print len(lines)
        while([s for s in lines[index].split(' ') if s!=''][0][0].isdigit()==False):
            #print [s for s in lines[index].split(' ') if s!='']
            #print index
            index+=1 
        while(lines[index].split(' ')[0]!=">>"):
            args = lines[index].split(' ')
            args = [s for s in args if s!='']
            if(args[0][0].isdigit()):
                new_record = HMMRecord(description = args[8],
                                                  e_value = float(args[0]),
                                                  score = float(args[1]),
                                                  bias = float(args[2]),
                                                  exp = float(args[6]),
                                                  num_domains = int(args[7]),
                                                  domains = [])
                #print>>debug_build,new_record
                self.hmm_records.append(new_record)
                
            index+=1
        self.num_records = len(self.hmm_records)
        #print self.num_records
    
    """
    Creates a file with the adapter sequences masked
    """
    def cleave_adapters(self,fasta_input_file,fasta_output_file):
        fin = open(fasta_input_file,'r')
        fout = open(fasta_output_file,'w')
        #debug = open('debug.txt','w')
        #debug1 = open('debug_data.txt','w')     
        lines = fin.readlines()
        hash = self.build_fasta_database(lines)
        for record in self.hmm_records:
            min_evalue = record.domains[0]
            min_dom = 0
            #Pick the more likely Domain (the smaller e-value)
            for i in range(0,len(record.domains)):
                if(record.domains[i].iEvalue<min_evalue):
                    min_evalue = record.domains[i].iEvalue
                    min_dom = i
            query_start = record.domains[min_dom].query_start
            query_end = record.domains[min_dom].query_end
            raw_seq = hash["> "+record.description+'\n']
            cleaved_seq = raw_seq[0:query_start-1]+\
                          raw_seq[query_start-1:query_end-1].lower()+\
                          raw_seq[query_end-1:]
            #Format all of the cleaved sequences into fasta format
            fout.write("> "+record.description+'\n')
            fout.write(cleaved_seq)
    """
    Creates a hashtable for the fasta input file
    """     
    def build_fasta_database(self,lines):
        
        database = dict([(lines[index],lines[index+1]) for index in range(0,len(lines),2)])
        return database
    
class HMMRecord(object):
    """
    Records information about each hit
    """
    def __init__(self,
                 description,
                 e_value,
                 score,
                 bias,
                 exp,
                 num_domains,
                 domains):

        """The constructor"""
        self.description = description
        self.e_value = e_value
        self.score = score
        self.bias = bias
        self.exp = exp
        self.num_domains = num_domains
        self.domains = domains
    def __str__(self):
        string = "%s \n Evalue: %f Score: %f Bias %f Exp %f"
        string = string%(self.description,self.e_value,self.score,self.bias,self.exp)
        for i in range(0,len(self.domains)):
            string+="\n\t==domain "+str(i+1)+"\n"+str(self.domains[i])
            i+=1
        return string
    def __del__(self):
        del self.domains

class Domain(object):
    """
    Records information about each potential alignment
    """
    def __init__(self,
                 score,
                 bias,
                 cEvalue,
                 iEvalue,
                 query_start,
                 query,
                 query_end,
                 sbjct_start,
                 sbjct,
                 sbjct_end):
        """The constructor"""
        self.score = score
        self.bias = bias
        self.cEvalue = cEvalue
        self.iEvalue = iEvalue
        self.query_start = query_start
        self.query = query
        self.query_end = query_end
        self.sbjct_start = sbjct_start
        self.sbjct = sbjct
        self.sbjct_end = sbjct_end
    def __str__(self):
        string =  "Score: "+str(self.score)+\
        "Bias: "+str(self.bias)+\
        "cEvalue: "+str(self.cEvalue)+\
        "iEvalue: "+str(self.iEvalue)+" \n"+\
            str(self.query_start)+"\t"+self.query+'\t'+str(self.query_end)+'\n'+\
            str(self.sbjct_start)+"\t"+self.sbjct+'\t'+str(self.sbjct_end)
        return string
                 
def read_chunks(file):
    fp = open(file)
    chunks_dict = {}
    for line in fp:
        line = line.rstrip()
        if line[0] == ">":
            line = line[1:]
            arr = re.split("\s+", line)
            chunks_dict[arr[0]] = arr[1]
    fp.close()
    return chunks_dict


def main():
    
    proc_id = sys.argv[2]
    chunks_dict = read_chunks(proc_id + "_chunks.fa")
    #H1 = HMMparser("hmmsearch.output")
    fp = open(sys.argv[1], 'r')
    fp_lines = fp.readlines()
    if len(fp_lines) == 36:            # In this case, there's no detected TE in the hmmsearch result file
        fa = open(proc_id + "_hmmparse.out", 'w')
        fa.write("The number of passed seqs is: 0\n")
    else:    
        H1 = HMMparser(sys.argv[1])

        fa = open(proc_id + "_hmmparse.out", 'w')
        records = H1.hmm_records
        fa.write("The number of passed seqs is: " + str(H1.num_records)+'\n')
        for i in range(len(records)):
            fa.write(records[i].description+'\t'+str(records[i].e_value)+'\t'+str(records[i].score)+'\n')
            dm = records[i].domains
            min_iEvalue = 1.0
            min_index = 0
            for j in range(len(dm)):
                if dm[j].iEvalue < min_iEvalue:
                    min_iEvalue = dm[j].iEvalue
                    min_index = j
            fa.write(str(dm[min_index].query_start+int(chunks_dict[records[i].description])-1)+'\t'+str(dm[min_index].query_end+int(chunks_dict[records[i].description]))+'\n')

        ##        if dm[j].iEvalue < 0.001:   # pick up the significant domains
        ##            fa.write(str(dm[j].query_start+int(chunks_dict[records[i].description])-1)+'\t'+str(dm[j].query_end+int(chunks_dict[records[i].description]))+'\n')
        fa.close()

if __name__=="__main__":
    main()
