import sys
import subprocess
import re
import os
import fileinput
import argparse
from glob import glob

proc_id = str(os.getpid())

def JC(a):    # a is a string number. Returns a list of 16 str(numbers) as the substitution rates
    subst_list = []
    b = str(1-3*float(a))
    subst_list.extend([b, a, a, a, a, b, a, a, a, a, b, a, a, a, a, b])
    return subst_list

def KM(a, b): # a, b are numbers in string format which are two args for the Kimura method. Returns a list of length 16
    subst_list = []
    c = str(1-float(a)-2*float(b))
    subst_list.extend([c, a, b, b, a, c, b, b, b, b, c, a, b, b, a, c])
    return subst_list
    

def one_simu(file_name, numFamilies, rng_seed = None):
##    str_randd = ""                    # change the randDis list to a string, will be a string like 'EXPO 1.0 5.0' or 'UNIF'
##    for i in range(len(randDis)):
##        str_randd = str_randd + randDis[i] + " "
##    cmd = "python one_whole_simu.py " + proc_id + " " + file_name + " " + numFamilies + str_randd + (" " + str(rng_seed) if rng_seed else "")
    cmd = "python one_whole_simu.py " + proc_id + " " + file_name + " " + numFamilies + (" " + str(rng_seed) if rng_seed else "")
    #print "one_whole_sim cmd: ", cmd
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();
    #print out

def set_params(subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd):
    """ change sim_params.txt
    all of the arguments here are [boolean, [value1, value2,...]]"""
    myfile = fileinput.input("sim_params.txt", inplace = 1)
    for line in myfile:
        line = line.rstrip()
        arr = re.split("\s+", line)   # split the line to an array
        if arr[0] == "substitution" and subst_rates[0] == True:  # substit rates
            print line
            line = myfile.next()
            line = subst_rates[1][0] + ' ' + subst_rates[1][1] + ' ' + subst_rates[1][2] + ' ' + subst_rates[1][3]
            print line
            line = myfile.next()
            line = subst_rates[1][4] + ' ' + subst_rates[1][5] + ' ' + subst_rates[1][6] + ' ' + subst_rates[1][7]
            print line
            line = myfile.next()
            line = subst_rates[1][8] + ' ' + subst_rates[1][9] + ' ' + subst_rates[1][10] + ' ' + subst_rates[1][11]
            print line
            line = myfile.next()
            line = subst_rates[1][12] + ' ' + subst_rates[1][13] + ' ' + subst_rates[1][14] + ' ' + subst_rates[1][15]
        elif arr[0] == "transition" and trans_rates[0] == True:  #transition rates
            print line
            line = myfile.next()
            line = trans_rates[1][0] + ' ' + trans_rates[1][1] + ' ' + trans_rates[1][2]
            print line
            line = myfile.next()
            line = trans_rates[1][3] + ' ' + trans_rates[1][4] + ' ' + trans_rates[1][5]
            print line
            line = myfile.next()
            line = trans_rates[1][6] + ' ' + trans_rates[1][7] + ' ' + trans_rates[1][8]
        elif arr[0] == "ancestor" and arr[1] == "distribution" and ancestor_dist[0] == True:  #ances dists
            print line
            line = myfile.next()
            line = ancestor_dist[1][0] + ' ' + ancestor_dist[1][1] + ' ' + ancestor_dist[1][2] + ' ' + ancestor_dist[1][3]
        elif arr[0] == "background" and arr[1] == "distribution" and backg_dist[0] == True:   #backg dists
            print line
            line = myfile.next()
            line = backg_dist[1][0] + ' ' + backg_dist[1][1] + ' ' + backg_dist[1][2] + ' ' + backg_dist[1][3]
        elif arr[0] == "ancestor" and arr[1] == "length" and ancestor_length[0] == True: #ances_length
            print line
            line = myfile.next()
            line = ancestor_length[1]
        elif arr[0] == "number" and arr[1] == "modern" and num_descend[0] == True:         #number of descendents
            print line
            line = myfile.next()
            line = num_descend[1]
        elif arr[0] == "cut" and cut_off[0] == True:         # Probability of cut-off for each TE(descendent) in the genome
            print line
            line = myfile.next()
            line = cut_off[1]
        elif arr[0] == "number" and arr[1] == "TEs" and num_TEs[0] == True:         # Number of TEs in the genome
            print line
            line = myfile.next()
            line = num_TEs[1]
        elif arr[0] == "background" and arr[1] == "length" and backg_len[0] == True:         # Background length between TEs
            print line
            line = myfile.next()
            line = backg_len[1]
        if arr[0] == "distribution" and randd[0] == True:         # Distribution to generate random numbers for TE cut off
            print line
            line = myfile.next()
            str_randd = ""                    # change the random command list to a string
            for i in range(len(randd[1])):
                str_randd = str_randd + randd[1][i] + " "
            line = str_randd
        print line

    fileinput.close()


def main():
    subst_rates = [False, []]
    trans_rates = [False, []]
    ancestor_dist = [False, []]
    backg_dist = [False, []]
    ancestor_length = [False, None]
    num_descend = [False, None]
    cut_off = [False, None]
    num_TEs = [False, None]
    backg_len = [False, None]
    randd = [False, []]

    parser = argparse.ArgumentParser(description = "Change parameters in sim_params.txt")
    parser.add_argument('-s', '--subst_rates', action = 'store', dest = 'subst', nargs = '+', \
                        default = [], type = str, help = 'New substitution rate matrix (A, C, G, T)')
    parser.add_argument('-t', '--trans_rates', action = 'store', dest = 'trans', nargs = 9, \
                        default = [], type = str, help = 'New transition rate matrix (Match, Insert, Delete)')
    parser.add_argument('-ad', '--ancestor_dist', action = 'store', dest = 'ances', nargs = 4, \
                        default = [], type = str, help = 'New background distribution of bases for the ancestor')
    parser.add_argument('-bd', '--backg_dist', action = 'store', dest = 'backg', nargs = 4, \
                        default = [], type = str, help = 'New background distribution of bases for the genome')
    parser.add_argument('-al', '--ancestor_length', action = 'store', dest = 'anlen', default = None, \
                        type = str, help = 'New length of ancestor')
    parser.add_argument('-nd', '--num_descend', action = 'store', dest = 'numdes', default = None, \
                        type = str, help = 'New number of descendents in the mal file')
    parser.add_argument('-co', '--cut_off', action = 'store', dest = 'cut', default = None, \
                        type = str, help = 'Probability of cut-off for each TE(descendent) in the genome')
    parser.add_argument('-nt', '--num_TEs', action = 'store', dest = 'numTEs', default = None, \
                        type = str, help = 'Number of TEs in the genome')
    parser.add_argument('-bl', '--backg_len', action = 'store', dest = 'backgl', default = None, \
                        type = str, help = 'Background length between TEs')
    parser.add_argument('-ph', '--print_header', action = 'store_true', default = False, help = \
                        'Only print header in the output file without running the whole process')
    parser.add_argument('-fn', '--file_name', action = 'store', dest = 'fname', default = 'analyserResult.txt', \
                        type = str, help = 'Name of the output file')
    parser.add_argument('--seed', action = 'store', dest = 'rng_seed', default = None, \
                        type = int, help = 'RNG Seed')
    parser.add_argument('-nf', '--num_family', action = 'store', dest = 'numFamilies', default = '1', \
                        type = str, help = 'Number of different families that we want to generate with the simulator')
    parser.add_argument('-dl', '--delete_files', action = 'store_true', default = False, help = \
                        'Delete all the files that has been generated except the result file')
    parser.add_argument('-rd', '--randd', action = 'store', dest = 'randDis', nargs = '+', \
                        default = [], type = str, help = 'The way to generate random numbers for the TE cut off probability in simulator, using uniform(one arg) or expo(two args)')


    args = parser.parse_args()
    #print "wrapper.py: ", arg.subst

    subst_list = []
    if args.subst != []:
        #print "args.subst: "+str(args.subst)
        subst_rates[0] = True
        subst_list = eval(args.subst[0])(*args.subst[1:])
        #if args.subst[0] == "JC":                           # JC method
        #    subst_list = substJC(*args.subst[1:])
        #elif args.subst[0] == "KM":                      # Kimura method with 3 args
        #    subst_list = substKM(*args.subst[1:])
            
    if args.trans != []:
        trans_rates[0] = True
    if args.ances != []:
        ancestor_dist[0] = True
    if args.backg != []:
        backg_dist[0] = True
    if args.anlen != None:
        ancestor_length[0] = True
    if args.numdes != None:
        num_descend[0] = True
    if args.cut != None:
        cut_off[0] = True
    if args.numTEs != None:
        num_TEs[0] = True
    if args.backgl != None:
        backg_len[0] = True
    if args.randDis != []:
        randd[0] = True
    
    subst_rates[1] = subst_list
    trans_rates[1] = args.trans
    ancestor_dist[1] = args.ances
    backg_dist[1] = args.backg
    ancestor_length[1] = args.anlen
    num_descend[1] = args.numdes
    cut_off[1] = args.cut
    num_TEs[1] = args.numTEs
    backg_len[1] = args.backgl
    randd[1] = args.randDis
    #print subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd
        
    set_params(subst_rates, trans_rates, ancestor_dist, backg_dist, ancestor_length, num_descend, cut_off, num_TEs, backg_len, randd)

    if args.print_header:
        f = open(args.fname, 'w')
        f.write('TP'.ljust(10)+'FP'.ljust(10)+'TN'.ljust(10)+'FN'.ljust(10)+'sensitivity'.ljust(20)+'averageSquareTrim'.ljust(20)+ \
                'leftMedian'.ljust(20)+'rightMedian'.ljust(20)+'leftMedianPositive'.ljust(20)+'rightMedianPositive'.ljust(20)+ \
                'leftMedianNegative'.ljust(20)+'rightMedianNegative\n')
        f.close()
    else:
        one_simu(args.fname, args.numFamilies, args.rng_seed)

    if args.delete_files:
        for f in glob ('*%s*'%proc_id):
            os.remove(f)



if __name__=="__main__":
    main()
