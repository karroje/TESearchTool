#!/opt/local/bin/python
import re
import sys
import random
import os
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

        else:
            continue
        
    fileptr.close()
    return [subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd]

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


def generate_descendts(cut_ends, mal_output, fa_output, true_coords, proc_id, num_family):

    track_coord = 0    #to record the true coordinate in the genome for each family

    for k in range(int(num_family)):
        ancst = generate_ancestor()
        ancestor_str = "".join(ancst)
        descendts = []     #to store all the descendents
        candi_TEs = []     #100 extra descendents for creating chromosome, chrom contains 100 TEs
        true_coord_pairs = []   #store the true coordinate pairs of the TEs in the genome
        temp_true_coord = backg_len + track_coord     #to record the coordinate of current TE (left side)
        for j in range(num_descend + num_TEs):   # num_descend: for .mal file to generate .hmm; num_TEs: for .fa file
            temp_dec = []            #to store temp descendent
            temp_dec_coord = []      #to store temp descendent with coordinates

            temp_rand = random.random()
            if j >= num_descend and temp_rand < cut_off:        #eg. if cut_off = .2, then 20% prob for this TE to be cut off
                if randd[0] == 'UNIF':
                    start_posi = random.randint(0,ancestor_length-1)   #start position in ancst
                    end_posi = random.randint(0,ancestor_length-1)     #close interval on both sides
                    if start_posi > end_posi:   #swap the start and end if start > end
                        temp = end_posi
                        end_posi = start_posi
                        start_posi = temp
                elif randd[0] == 'EXPO':
                    cutoffExpo = cutoffFromExpo(float(randd[1]), float(randd[2]), ancestor_length)
                    start_posi = cutoffExpo[0]
                    end_posi = cutoffExpo[1]
            else:
                start_posi, end_posi = 0, ancestor_length-1
            temp_length = end_posi - start_posi + 1 # length of current descendent
            i = 0
            state = 0      # state = 0 if match, 1 if insert, 2 if delete
            while i < temp_length:
                if i == 0:     # the first letter should only be match
                    temp_dec.append(BASES[pickFromDistribution(subst_rates[base2index[ancst[start_posi]]])])
                    i += 1
                    state = 0
                else:   # if this is not the first letter, then match/insert/delete may happen
                    if state == 0:    # if the state of the previous letter is match, then:
                        if pickFromDistribution(trans_rates[0]) == 0:    # from match to match
                            ancst_base_index = base2index[ancst[start_posi + i]]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                            i += 1
                            state = 0
                        elif pickFromDistribution(trans_rates[0]) == 1:   # from match to insert
                            temp_dec.append(BASES[pickFromDistribution(backg_dist)].lower())
                            state = 1
                        elif pickFromDistribution(trans_rates[0]) == 2:   # from match to delete      
                            temp_dec.append("-")
                            i += 1
                            state = 2
                    elif state == 1:  # if the state of the previous letter is insert, then:
                        if pickFromDistribution(trans_rates[1]) == 0:     # from insert to match
                            ancst_base_index = base2index[ancst[start_posi + i]]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])
                            i += 1
                            state = 0
                        elif pickFromDistribution(trans_rates[1]) == 1:   # from insert to insert
                            temp_dec.append(BASES[pickFromDistribution(backg_dist)].lower())
                            state = 1
                    elif state == 2:  # if the state of the previous letter is delete, then:
                        if pickFromDistribution(trans_rates[2]) == 0:     # from delete to match
                            ancst_base_index = base2index[ancst[start_posi + i]]
                            temp_dec.append(BASES[pickFromDistribution(subst_rates[ancst_base_index])])                        
                            i += 1
                            state = 2
                        elif pickFromDistribution(trans_rates[2]) == 2:   # from delete to delete
                            temp_dec.append("-")
                            i += 1
                            state = 2
                    else:
                        print "Wrong state!"
            temp_dec = "".join(temp_dec)
            temp_dec_coord.append(temp_dec)
            temp_dec_coord.append(start_posi)
            temp_dec_coord.append(end_posi + 1)   # half open interval

            if j < num_descend:
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
        for i in range(num_TEs):
        #for i in range(len(true_coord_pairs)):
            fileptr_w4.write(str(true_coord_pairs[i][0]) + '\t' + str(true_coord_pairs[i][1]) + '\n')
        fileptr_w4.close()

        track_coord = temp_true_coord - backg_len
        
     

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

    subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd = read_params(sys.argv[marker])
    #generate_descendts(generate_ancestor(), cut_ends, sys.argv[marker+1], sys.argv[marker+2], sys.argv[marker+3], sys.argv[marker+4])
    generate_descendts(cut_ends, sys.argv[marker+1], sys.argv[marker+2], sys.argv[marker+3], sys.argv[marker+4], sys.argv[marker+5])

    
