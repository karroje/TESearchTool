#!/opt/local/bin/python
import re
import sys
import random

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
                ancestor_dist.append(temp_cdf[i])
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
            ancestor_length = arr[0]

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
    return [subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, gen_length, start_in_gen]

def pickFromDistribution(dist):
    b = random.random()
    i = 0;
    while (i < len(dist)-1 and b > dist[i]):
        i += 1
    return i

def generate_ancestor():
    ancestor = []
    for i in range(ancestor_length):
        ancestor.append(BASES[pickFromDistribution(ancestor_dist)])
    return ancestor                 #ancestor is a list here(not a str)


def generate_descendts(ancst, cut_ends, mal_output, fa_output):
    ancestor_str = "".join(ancst)
    descendts = []     #to store all the descendents
    for j in range(num_descend + 1):   # this "+1" is to generate an extra descendent for .fa file
        temp_dec = []            #to store temp descendent
        temp_dec_coord = []      #to store temp descendent with coordinates

        if cut_ends:
            start_posi = random.randint(0,ancestor_length-1)   #start position in ancst
            end_posi = random.randint(0,ancestor_length-1)     #close interval on both sides
            if start_posi > end_posi:   #swap the start and end if start > end
                temp = end_posi
                end_posi = start_posi
                start_posi = temp
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

        if j != num_descend:
            descendts.append(temp_dec_coord)
        else:                                 # generate an extra descendent for .fa file
            temp_for_fa_file = temp_dec_coord    

    descend_with_backg = []       # descendent with background
    for i in range(start_in_gen):
        descend_with_backg.append(BASES[pickFromDistribution(backg_dist)])
    descend_with_backg.append(temp_for_fa_file[0])
    for i in range(gen_length - start_in_gen):
        descend_with_backg.append(BASES[pickFromDistribution(backg_dist)])
    descend_with_backg = "".join(descend_with_backg)
    descend_with_backg = descend_with_backg.upper()
    fileptr_w2 = open(fa_output, "w")
    fileptr_w2.write(">" + str(start_in_gen) + "\t" + str(temp_for_fa_file[1]) + "\t" + str(temp_for_fa_file[2]) + "\t" + str(temp_for_fa_file[0]) + "\n" + "".join([x for x in descend_with_backg if x!='-']))
        
    fileptr_w1 = open(mal_output, "w")
    fileptr_w1.write("ancestor" + "\t" + ancestor_str + "\n")
    for i in range(len(descendts)):
        fileptr_w1.write("chr?" + "\t" + "?" + "\t" + "?" + "\t" + str(descendts[i][1]) + "\t" + str(descendts[i][2]) + "\t" + str(descendts[i][0]) + "\n")
    fileptr_w1.close()

    fileptr_w3 = open("SimuRepBase.fa", "w")
    fileptr_w3.write(">SimuRepBase" + '\n' + ancestor_str.lower())
    fileptr_w3.close()
        
     

if __name__=="__main__":
    cut_ends = False
    seed = 0

    marker = 1
    while marker < len(sys.argv) and sys.argv[marker][0] == "-":
        s = sys.argv[marker]
        if s == "-c" or s == "-cut":
            cut_ends = True
        elif s == "-s":
            marker = marker + 1;
            rng = int(sys.argv[marker])
        marker = marker + 1
    assert(marker == len(sys.argv)-3)

    if seed != 0:
        random.seed(seed)

    subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, gen_length, start_in_gen = read_params(sys.argv[marker])
    generate_descendts(generate_ancestor(), cut_ends, sys.argv[marker+1], sys.argv[marker+2])

    
