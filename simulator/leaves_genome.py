import sys
import re
import random
from tree import *

BASES = ['A', 'C', 'G', 'T']
base2index = {'A':0, 'C':1, 'G':2, 'T':3, 'a':0, 'c':1, 'g':2, 't':3}

#read parameters, matrices from sim_params.txt
def read_params(file):
    """Reed in the simulation parameters from the speficied file."""
    subst_rates = []     
    trans_rates = []
    ances_dist = []
    backg_dist = []
    ances_length = 0
    num_descend = 0
    gen_length = 0
    start_in_gen = 0
    
    fileptr = open(file)
    for line in fileptr:
        if line[0] == '#':    #skip comments
            next;    
        line = line.rstrip()
        arr = re.split("\s+", line)   # split the line to an array

        if arr[0] == "substitution":  # substit rates
            line = fileptr.next()     # is this correct??
            for i in range(4):        # 4 lines in this matrix
                arr = re.split("\s+", line)
                temp_cdf = [0.0,0.0,0.0,0.0]
                for i in range(4):    # 4 numbers in each line
                    arr[i] = float(arr[i])   # change the probs to floats
                    temp_cdf[i] = temp_cdf[i-1] + arr[i]   #calculate cdf
                arr = temp_cdf
                subst_rates.append(arr)
                line = fileptr.next()

        elif arr[0] == "transition":  #transition rates
            line = fileptr.next()     
            for i in range(3):        # 3 lines in this matrix
                arr = re.split("\s+", line)
                temp_cdf = [0.0,0.0,0.0]
                for i in range(3):    # 3 numbers in each line
                    arr[i] = float(arr[i])   # change the probs to floats
                    temp_cdf[i] = temp_cdf[i-1] + arr[i]   #calculate cdf
                arr = temp_cdf
                trans_rates.append(arr)
                line = fileptr.next()

        elif arr[0] == "ancestor" and arr[1] == "distribution":  #ances dists
            line = fileptr.next()
            arr = re.split("\s+", line)
            temp_cdf = [0.0,0.0,0.0,0.0]
            for i in range(4):    
                arr[i] = float(arr[i])   # change the probs to floats
                temp_cdf[i] = temp_cdf[i-1] + arr[i]   #calculate cdf
                ances_dist.append(temp_cdf[i])
            line = fileptr.next()    

        elif arr[0] == "background":          #backg dists
            line = fileptr.next()
            arr = re.split("\s+", line)
            temp_cdf = [0.0,0.0,0.0,0.0]
            for i in range(4):    
                arr[i] = float(arr[i])   # change the probs to floats
                temp_cdf[i] = temp_cdf[i-1] + arr[i]   #calculate cdf
                backg_dist.append(temp_cdf[i])
            line = fileptr.next()

        elif arr[0] == "ancestor" and arr[1] == "length": #ances_length
            line = fileptr.next()
            arr = re.split("\s+", line)
            arr[0] = int(arr[0])
            ances_length = arr[0]

        elif arr[0] == "number":         #number of descendents
            #print "X: ", num_descend
            line = fileptr.next()
            arr = re.split("\s+", line)
            #print arr
            arr[0] = int(arr[0])
            num_descend = arr[0]
            #print num_descend

        elif arr[0] == "genome":         #genome length
            line = fileptr.next()
            arr = re.split("\s+", line)
            arr[0] = int(arr[0])
            gen_length = arr[0]

        elif arr[0] == "repeat":         #repeat position
            line = fileptr.next()
            arr = re.split("\s+", line)
            arr[0] = int(arr[0])
            start_in_gen = arr[0]

        else:
            continue
        
    fileptr.close()
    return [subst_rates, trans_rates, ances_dist, backg_dist, ances_length, num_descend, gen_length, start_in_gen]


def pickFromDistribution(dist):
    b = random.random()
    i = 0;
    while (i < len(dist)-1 and b > dist[i]):
        i += 1
    return i

def genrt_ances():
    ancestor = []
    for i in range(ances_length):
        ancestor.append(BASES[pickFromDistribution(ances_dist)])
    return ancestor                 #ancestor is a list here(not a str)
    

def genrt_descendts(ancst, n):  # n is the number of times of mutation from the ancestor
    ances_str = "".join(ancst)
    for j in range(n):
        temp_dec = []            #to store temp descendent
        temp_dec_coord = []      #to store temp descendent with coordinates
        start_posi = 0           #full length
        end_posi = len(ances_str) - 1
        temp_length = len(ances_str)
        i = 0
        state = 0      # state = 0 if match, 1 if insert, 2 if delete
        while i < temp_length:
            if i == 0:     # the first letter should only be match
                temp_dec.append(BASES[pickFromDistribution(subst_rates[base2index[ances_str[start_posi]]])])
                i += 1
                state = 0
            else:   # if this is not the first letter, then match/insert/delete may happen
                if ances_str[start_posi + i] == '-':    # if it's dash, just copy it and increase index
                    temp_dec.append("-")
                    i += 1
                else:
                    if state == 0:    # if the state of the previous letter is match, then:
                        new_state = pickFromDistribution(trans_rates[0])
                        if new_state == 0:    # from match to match
                            ancst_base_index = base2index[ances_str[start_posi + i]]
                            if ances_str[start_posi + i].islower():
                                temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])].lower())
                            else:
                                temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                            i += 1
                        elif new_state == 1:   # from match to insert
                            temp_dec.append(BASES[pickFromDistribution(backg_dist)].lower())
                        else: # new_state == 2   from match to delete      
                            temp_dec.append("-")
                            i += 1
                        state = new_state
                        
                    elif state == 1:  # if the state of the previous letter is insert, then:
                        if pickFromDistribution(trans_rates[1]) == 0:     # from insert to match
                            ancst_base_index = base2index[ances_str[start_posi + i]]
                            if ances_str[start_posi + i].islower():
                                temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])].lower())
                            else:
                                temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                            i += 1
                            state = 0
                        else:   # from insert to insert
                            temp_dec.append(BASES[pickFromDistribution(backg_dist)].lower())
                            state = 1
                    elif state == 2:  # if the state of the previous letter is delete, then:
                        if pickFromDistribution(trans_rates[2]) == 0:     # from delete to match
                            ancst_base_index = base2index[ances_str[start_posi + i]]
                            if ances_str[start_posi + i].islower():
                                temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])].lower())
                            else:
                                temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])]) 
                            i += 1
                            state = 2
                        else:   # from delete to delete
                            temp_dec.append("-")
                            i += 1
                            state = 2
                    else:
                        print "Wrong state!"

        temp_dec = "".join(temp_dec)
        #print temp_dec
        ances_str = temp_dec
    return temp_dec


def genrt_node_seqs(t):      # to generate sequences for each node of tree t, need the root sequence assigned
    if t.isLeaf() == False:      # while it has children
        for i in range(t.numChildren()):
            t.child(i)[0].sequence = genrt_descendts(t.sequence, int(t.child(i)[1]))
            genrt_node_seqs(t.child(i)[0])
        return t.sequence
    else:
        return t.sequence
        
def genrt_whole_seq(seq_list, back_len):   # to generate the whole seq with background
    for i in range(len(seq_list)):         # to reduce dashes and make them all upper cases
        seq_list[i] = seq_list[i].upper()
        seq_list[i] = re.split('-', seq_list[i])
        seq_list[i] = "".join(seq_list[i])
    whole_seq = []             # sequence line
    whole_coords = []          # coordinates of each descendent
    for i in range(len(seq_list)):      # to calculate whole_coords
        seq_with_coord = []
        if i == 0:
            temp_start = back_len   # back_len is the background length between any two neighbor descendents
        else:
            temp_start = whole_coords[i-1][2] + back_len    # whole_coords[i-1][2] is the end coord of the previous seq
        temp_end = temp_start + len(seq_list[i])
        seq_with_coord.append(seq_list[i])
        seq_with_coord.append(temp_start)
        seq_with_coord.append(temp_end)
        whole_coords.append(seq_with_coord)
    for j in range(len(seq_list)):      # to make whole_seq with background
        for k in range(back_len):
            whole_seq.append(BASES[pickFromDistribution(backg_dist)])
        whole_seq.append(seq_list[j])
        if j == len(seq_list)-1:
            for k in range(back_len):
                whole_seq.append(BASES[pickFromDistribution(backg_dist)])
    whole_seq = "".join(whole_seq)
    return [whole_seq, whole_coords]

def gener_output_file(whole_thing, out):   # the input format is like the result of def genrt_whole_seq
    fileptr_w1 = open(out + ".fa", "w")
    fileptr_w1.write(whole_thing[0] + "\n" + ">")   # the whole seq
    for i in range(len(whole_thing[1])):      # the coords of each descendent in the seq
        fileptr_w1.write(str(whole_thing[1][i][1]) + "\t" + str(whole_thing[1][i][2]) + "\t")
    fileptr_w1.close()
        
        


new_str = "(A:10)R"
#new_str = "((D:5,E:2)C:3,B:4)A"
t = parseNewick(new_str)
##print "countLeaves:", t.countLeaves()
##print "totalWeight:", t.totalWeight()
##print "leafList: ", t.leafList()
##print "nodeDepths: ", t.nodeDepths()
##print "printNewick: ", t.printNewick()
##
##str2 = splitTopLevel(new_str)
##print str2

subst_rates, trans_rates, ances_dist, backg_dist, ances_length, num_descend, gen_length, start_in_gen = read_params("sim_params.txt")
node_ances = "".join(genrt_ances())    # create the sequence of the root
t.set_seq(node_ances)                  # set the sequence of the root
#print t.sequence
genrt_node_seqs(t)
sequence_list = t.leafSeqs()
#print genrt_whole_seq(sequence_list, 100)
whole_thing = genrt_whole_seq(sequence_list, 100)
gener_output_file(whole_thing, "output")
