#!/opt/local/bin/python
import re
import sys
import random
import os
from tree import *
#from multi_simu import proc_id

"""Read in the sim_params.txt file and generate simulated data: 
1. _sim_out.mal : consists of the ancestor and descendents; will be converted to .align2 file for the TEprofiler
2. _SimuRepBase.fa : containing the ancestor as the family sequence; will be used by TEprofiler and chunks.cpp
3. _sim_out.fa : the genome which contains multiple TEs; will be used by hunter and chunks.cpp
4. _true_coords.txt: to record the TEs' true coordinates; will be used for the result_analyser.py
"""

BASES = ['A', 'C', 'G', 'T']
base2index = {'A':0, 'C':1, 'G':2, 'T':3}

def read_params(file):
    """Reed in the simulation parameters from the speficied file."""
    subst_rates = []     
    trans_rates = []
    ancestor_dist = []
    backg_dist = []
    ancestor_length = 0
    num_descend = 0
##    gen_length = 0
##    start_in_gen = 0
    cut_off = 0.0
    num_TEs = 0
    backg_len = 0
    randd = []
    tree_str = ""
    
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
                ancestor_dist.append(temp_cdf[i])
            line = fileptr.next()    

        elif arr[0] == "background" and arr[1] == "distribution":          #backg dists
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
            ancestor_length = arr[0]

        elif arr[0] == "number" and arr[1] == "modern":         #number of descendents
            #print "X: ", num_descend
            line = fileptr.next()
            arr = re.split("\s+", line)
            #print arr
            arr[0] = int(arr[0])
            num_descend = arr[0]
            #print num_descend

        elif arr[0] == "cut":            #Probability of cut-off for each TE(descendent) in the genome
            line = fileptr.next()
            arr = re.split("\s+", line)
            arr[0] = float(arr[0])
            cut_off = arr[0]

        elif arr[0] == "number" and arr[1] == "TEs":    # Number of TEs in the genome
            line = fileptr.next()
            arr = re.split("\s+", line)
            arr[0] = int(arr[0])
            num_TEs = arr[0]

        elif arr[0] == "background" and arr[1] == "length":   # Background length between TEs
            line = fileptr.next()
            arr = re.split("\s+", line)
            arr[0] = int(arr[0])
            backg_len = arr[0]

        elif arr[0] == "distribution":   # Distribution to generate random numbers for TE cut off
            line = fileptr.next()
            arr = re.split("\s+", line)
            randd = arr

        elif arr[0] == "tree":   # Tree discription of the TE family (we will use either this option or "number TEs in the genome"
            line = fileptr.next()
            tree_str = line.rstrip()

        else:
            continue
        
    fileptr.close()
    return [subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd, tree_str]

def pickFromDistribution(dist):
    b = random.random()
    i = 0;
    while (i < len(dist)-1 and b > dist[i]):
        i += 1
    return i

def cutoffFromExpo(mean1, mean2, lenTE):   # to generate the start and end position in one TE using exponential dist
    start = end = -1
    while(not(start >=0 and end >=0 and start + end < lenTE)):
        start = int(round(random.expovariate(mean1)))
        end = int(round(random.expovariate(mean2)))
    start_coord = start
    end_coord = lenTE - end - 1
    return [start_coord, end_coord]

def generate_ancestor():
    ancestor = []
    for i in range(ancestor_length):
        ancestor.append(BASES[pickFromDistribution(ancestor_dist)])
    return ancestor                 #ancestor is a list here(not a str)

def gener_single_des(ancst, start_posi, end_posi):
    #print subst_rates
    temp_dec = []            #to store temp descendent
    temp_length = end_posi - start_posi + 1 # length of current descendent
    i = 0
    state = 0      # state = 0 if match, 1 if insert, 2 if delete
##    print "ancst_length: ", len(ancst)
##    print "templength: ", temp_length
##    print "tempacestor: ", ancst
##    print "temp_start: ", start_posi
##    print "temp_end: ", end_posi
    while i < temp_length:
        if i == 0:     # the first letter should only be match
            temp_dec.append(BASES[pickFromDistribution(subst_rates[base2index[ancst[start_posi]]])])
            i += 1
            state = 0
        else:   # if this is not the first letter, then match/insert/delete may happen
            if ancst[start_posi + i] == '-':   # if the current charater is dash, just copy it and define it as match state
                temp_dec.append('-')
                i += 1
                state = 0
            else:                              # if the current character is not dash, then:
                if state == 0:    # if the state of the previous letter is match, then:
                    if pickFromDistribution(trans_rates[0]) == 0:    # from match to match
                        if ancst[start_posi + i].isupper():          # if it's upper case
                            ancst_base_index = base2index[ancst[start_posi + i]]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                        else:                                        # if it's lower case (comes from insertions of previous iterations)
                            ancst_base_index = base2index[ancst[start_posi + i].upper()]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])].lower())
                        i += 1
                        state = 0
                    elif pickFromDistribution(trans_rates[0]) == 1:   # from match to insert
                        temp_dec.append(BASES[pickFromDistribution(backg_dist)].lower())
                        state = 1
                    elif pickFromDistribution(trans_rates[0]) == 2:   # from match to delete
                        if ancst[start_posi + i].isupper():           # if it's upper case, add dash; if not, do not add any char, just increase index
                            temp_dec.append("-")
                        i += 1
                        state = 2
                elif state == 1:  # if the state of the previous letter is insert, then:
                    if pickFromDistribution(trans_rates[1]) == 0:     # from insert to match
                        if ancst[start_posi + i].isupper():           # if it's upper case
                            ancst_base_index = base2index[ancst[start_posi + i]]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                        else:                                         # if it's lower case (comes from insertions of previous iterations)
                            ancst_base_index = base2index[ancst[start_posi + i].upper()]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])].lower())
                        i += 1
                        state = 0
                    elif pickFromDistribution(trans_rates[1]) == 1:   # from insert to insert
                        temp_dec.append(BASES[pickFromDistribution(backg_dist)].lower())
                        state = 1
                elif state == 2:  # if the state of the previous letter is delete, then:
                    if pickFromDistribution(trans_rates[2]) == 0:     # from delete to match
                        #print i, start_posi+i
                        if ancst[start_posi + i].isupper():           # if it's upper case
                            ancst_base_index = base2index[ancst[start_posi + i]]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                        else:                                         # if it's lower case (comes from insertions of previous iterations)
                            ancst_base_index = base2index[ancst[start_posi + i].upper()]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])].lower())                        
                        i += 1
                        state = 0
                    elif pickFromDistribution(trans_rates[2]) == 2:   # from delete to delete
                        temp_dec.append("-")
                        i += 1
                        state = 2
                else:
                    print "Wrong state!"
    temp_dec = "".join(temp_dec)
    return temp_dec


"""Change the subst_rates according to the alpha on the edge of the tree
param_list looks like: [100, 0.01] or [100, 0.01, 0.02], the number of params determines
which method to use, 2 for JC, 3 for Kimura
"""
def changeSubst(param_list):
    #print param_list
    new_subst_pdf = []
    if len(param_list)==2:   # JC
        a = float(param_list[1])
        b = 1.0-3*a
        new_subst_pdf.append([b, a, a, a])
        new_subst_pdf.append([a, b, a, a])
        new_subst_pdf.append([a, a, b, a])
        new_subst_pdf.append([a, a, a, b])
    elif len(param_list)==3:    #Kimura
        a = float(param_list[1])
        b = float(param_list[2])
        c = 1.0-a-2*b
        new_subst_pdf.append([c, a, b, b])
        new_subst_pdf.append([a, c, b, b])
        new_subst_pdf.append([b, b, c, a])
        new_subst_pdf.append([b, b, a, c])
    #print new_subst_pdf
    new_subst_cdf = []
    for i in range(4):        # 4 lines in this matrix
        #print i
        arr = new_subst_pdf[i]
        temp_cdf = [0.0,0.0,0.0,0.0]
        for i in range(4):    # 4 numbers in each line
            temp_cdf[i] = temp_cdf[i-1] + arr[i]   #calculate cdf
        arr = temp_cdf
        new_subst_cdf.append(arr)
    return new_subst_cdf


"""Generates sequences for each node of tree t, need the root sequence assigned"""
def genrt_node_seqs(t):          
    if t.isLeaf() == False:                      # while it has children
        for i in range(t.numChildren()):
            global subst_rates
            subst_rates = changeSubst(re.split('\|', t.child(i)[1]))    # change the subst_rates
            #t_seq = "".join(re.split('-', t.sequence)).upper()
            #t.child(i)[0].sequence = gener_single_des(t_seq, 0, len(t_seq)-1)      # the seqs for interior nodes have no cut off
            t.child(i)[0].sequence = gener_single_des(t.sequence, 0, len(t.sequence)-1)
            #t.child(i)[0].base_seq = t_seq
            t.child(i)[0].base_seq = t.sequence
            genrt_node_seqs(t.child(i)[0])
        return t.sequence
    else:
        return t.sequence


"""This method is for the simulator, returns a list of the leaf nodes.
Each element in the list corresponds to a leaf, which contains:
1. a sequence that is the same as its parent interior node
2. params list such as [100, 0.01]
"""
def listForSimu(t):
    list_for_simu = []
    listForSimuHelper(t, list_for_simu)
    return list_for_simu

def listForSimuHelper(t, list_for_simu):
    if t.isLeaf():
        return list_for_simu
    else:
        for i in range(t.numChildren()):
            if t.child(i)[0].isLeaf():              # if this child tree is a leaf, add the subList to the list_for_simu
                list_for_simu.append([t.child(i)[0].base_seq, re.split('\|', t.child(i)[1])])     # child(i)[1] is the params list for child(i) such as: [100, 0.01]
            else:
                listForSimuHelper(t.child(i)[0], list_for_simu)
        return list_for_simu


def generate_descendts(cut_ends, mal_output, fa_output, true_coords, proc_id, num_family):

    track_coord = 0    #to record the true coordinate in the genome for each family

    for k in range(int(num_family)):
        ancst = generate_ancestor()
        ancestor_str = "".join(ancst)

        t.set_seq(ancestor_str)         # t is a global variable
        '''do the recursion with the tree to get a list for the leaves'''
        genrt_node_seqs(t)
        list_for_simu = listForSimu(t)  # looks like [['AAATT', ['100', '0.01']], ['TTCCA', ['200', '0.02']]]
        
        descendts = []     #to store all the descendents for the .mal file to generate .hmm;
        candi_TEs = []     #extra descendents for creating chromosome, chrom contains TEs, .fa file
        true_coord_pairs = []   #store the true coordinate pairs of the TEs in the genome
        temp_true_coord = backg_len + track_coord     #to record the coordinate of current TE (left side)
        
        #for j in range(num_descend + num_TEs):   # num_descend: for .mal file to generate .hmm; num_TEs: for .fa file
            
        for h in range(len(list_for_simu)):
            temp_ancst = list_for_simu[h][0]                # ancestor sequence for current iteration
            print "len_temp_ancst: ", len(temp_ancst)
            if len(temp_ancst)==501:
                print temp_ancst[0], temp_ancst[500]
            temp_ancestor_length = len(temp_ancst)
            num_for_this_leaf = int(list_for_simu[h][1][0])      # number of descendents needed to generate for this leaf
            global subst_rates
            subst_rates = changeSubst(list_for_simu[h][1])  # change the substitution rate matrix based on the list of parameters for the current leaf
            
            for j in range(num_for_this_leaf * 2):          # two parts of TEs, first half for .mal, second half for .fa
            
                temp_dec_coord = []      #to store temp descendent with coordinates

                temp_rand = random.random()
                if temp_rand < cut_off:        #eg. if cut_off = .2, then 20% prob for this TE to be cut off
                    if randd[0] == 'UNIF':
                        start_posi = random.randint(0,temp_ancestor_length-1)   #start position in ancst
                        end_posi = random.randint(0,temp_ancestor_length-1)     #close interval on both sides
                        if start_posi > end_posi:   #swap the start and end if start > end
                            temp = end_posi
                            end_posi = start_posi
                            start_posi = temp
                    elif randd[0] == 'EXPO':
                        cutoffExpo = cutoffFromExpo(float(randd[1]), float(randd[2]), temp_ancestor_length)
                        start_posi = cutoffExpo[0]
                        end_posi = cutoffExpo[1]
                else:
                    start_posi, end_posi = 0, temp_ancestor_length-1
                
                temp_dec = gener_single_des(temp_ancst, start_posi, end_posi)
                temp_dec_coord.append(temp_dec)
                temp_dec_coord.append(start_posi)
                temp_dec_coord.append(end_posi + 1)   # half open interval

                if j < num_for_this_leaf:
                    descendts.append(temp_dec_coord)
                else:                                 # generate TEs for .fa file
                    temp_for_fa_file = temp_dec_coord    
                    descend_with_backg = []       # descendent with background
                    for i in range(backg_len):
                        descend_with_backg.append(BASES[pickFromDistribution(backg_dist)])
                    descend_with_backg.append(temp_for_fa_file[0])
                    descend_with_backg = "".join(descend_with_backg)
                    descend_with_backg = descend_with_backg.upper()
                    candi_TEs.append(descend_with_backg)
                    
                    true_coord_pairs.append([temp_true_coord, temp_true_coord+len(temp_dec)-temp_dec.count('-')])
                    temp_true_coord = temp_true_coord + len(temp_dec) - temp_dec.count('-') + backg_len

        candi_TEs = "".join(candi_TEs)    
        fileptr_w2 = open(proc_id + "_" + fa_output, "a")
        if k == 0:
            fileptr_w2.write(">Simulated genome" + "\n" + "".join([x for x in candi_TEs if x!='-']))
        else:
            fileptr_w2.write("".join([x for x in candi_TEs if x!='-']))
            
        fileptr_w1 = open(str(k) + "_" + proc_id + "_" + mal_output, "w")
        fileptr_w1.write("ancestor" + "\t" + ancestor_str + "\n")
        for i in range(len(descendts)):
            fileptr_w1.write("chr?" + "\t" + "?" + "\t" + "?" + "\t" + str(descendts[i][1]) + "\t" + str(descendts[i][2]) + "\t" + str(descendts[i][0]) + "\n")
        fileptr_w1.close()

        fileptr_w3 = open(str(k) + "_" + proc_id + "_SimuRepBase.fa", "w")
        #fileptr_w3.write(">" + proc_id + "_SimuRepBase" + '\n' + ancestor_str.lower())
        fileptr_w3.write(">" + str(k) + "_" + proc_id + "_SimuRepBase" + '\n' + ancestor_str.lower())
        fileptr_w3.close()

        fileptr_w4 = open(str(k) + "_" + proc_id + "_" + true_coords, "w")
        print "num_TEs: ", num_TEs, "true_coord_pairs length: ", str(len(true_coord_pairs))
        #for i in range(num_TEs):
        for i in range(len(true_coord_pairs)):
            fileptr_w4.write(str(true_coord_pairs[i][0]) + '\t' + str(true_coord_pairs[i][1]) + '\n')
        fileptr_w4.close()

        track_coord = temp_true_coord - backg_len
        
     

##subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd = read_params("sim_params.txt")
##new_str = "((D:10|0.06,E:10|0.05)C:1|0.06,B:10|0.07)A"
##t = parseNewick(new_str)
##ancst = generate_ancestor()
##ancestor_str = "".join(ancst)
##t.set_seq(ancestor_str)
##print "ances:", t.sequence
##genrt_node_seqs(t)
##print "leafSeqs:", t.leafSeqs()
##print "leafBaseSeqs:", t.leafBaseSeqs()
##print "listForSimu:", listForSimu(t)

##print re.split('\|', t.child(1)[1])
##print "printNewick: ", t.printNewick()


if __name__=="__main__":
    cut_ends = False
    seed = 0

    marker = 1
    while marker < len(sys.argv) and sys.argv[marker][0] == "-":
        s = sys.argv[marker]
        if s == "-c" or s == "-cut":
            cut_ends = True
        elif s == "--seed":
            marker = marker + 1;
            seed = int(sys.argv[marker])
        marker = marker + 1
    assert(marker == len(sys.argv)-6)

    if seed != 0:
        random.seed(seed)

    subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd, tree_str = read_params(sys.argv[marker])
    #new_str = "((D:10|0.01,E:10|0.02)C:1|0.01,B:10|0.01)A"
    t = parseNewick(tree_str)
    generate_descendts(cut_ends, sys.argv[marker+1], sys.argv[marker+2], sys.argv[marker+3], sys.argv[marker+4], sys.argv[marker+5])

    
