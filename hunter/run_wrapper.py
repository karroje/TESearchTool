import sys
import subprocess
import re
import os
import fileinput
import argparse    

def add_arg_col(fname, arg_name_list, arg_val_list, numFamilies):
    '''add one colume on the left of the original result file, including the name and
    value. Just add the value to the last numFamilies lines if the name is already in the header'''
    lines = open(fname, 'r').readlines()
    '''For the first line'''
    line0_arr = re.split('\s+', lines[0].rstrip())
    if line0_arr[0]!=arg_name_list[0]:
        str_name_list = ""
        for i in range(len(arg_name_list)):
            str_name_list = str_name_list + arg_name_list[i].ljust(15)
        lines[0] = str_name_list + lines[0]
    '''For the last numFamilies lines'''
    str_val_list = ""
    for i in range(len(arg_val_list)):
        str_val_list = str_val_list + arg_val_list[i].ljust(15)
    for i in range(int(numFamilies)):
        lines[-1-i] = str_val_list + lines[-1-i]
    open(fname, 'w').writelines(lines)
    

def main():

    parser = argparse.ArgumentParser(description = "Set up the input arguments and invoke \"wrapper.py\"")
    
##    parser.add_argument('-s', '--subst_rates', action = 'store', dest = 'subst', default = None, \
##                        type = str, help = 'This subst is m, which defines the substitution matrix, first row is 1-3m, m, m, m')
##    parser.add_argument('-sj', '--subst_jc', action = 'store', dest = 'substjc', default = None, \
##                        type = str, help = 'New substitution rate matrix (A, C, G, T), using JC method with one parameters a. first row: 1-3a, a, a, a')
##    parser.add_argument('-sk', '--subst_km', action = 'store', dest = 'substkm', nargs = 2, \
##                        default = [], type = str, help = 'New substitution rate matrix (A, C, G, T), using Kimura method with two parameters a and b')
    parser.add_argument('-s', '--subst', action = 'store', dest = 'subst', nargs = '+', \
                        default = [], type = str, help = 'New substitution rate matrix (A, C, G, T), using JC(two args) or Kimura(three args)')

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
    #print "run_wrapper: ", args.fname, args.rng_seed

    if args.print_header:                # cannot be used with any other option
        cmd = "python wrapper.py -ph -fn " + args.fname
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();
    #elif args.subst != None:
    elif args.subst != []:

        str_subst = ""                    # change the list to a string
        for i in range(len(args.subst)):
            str_subst = str_subst + args.subst[i] + " "

        if args.randDis != []:
            str_randd = ""                    # change the random command list to a string
            for i in range(len(args.randDis)):
                str_randd = str_randd + args.randDis[i] + " "
            
        cmd = "python wrapper.py -s " + str_subst + ((" -rd " + str_randd) if args.randDis != [] else "") + " -fn "+args.fname + " -nf " +args.numFamilies + \
              (" -dl " if args.delete_files else "") + ((" --seed " + str(args.rng_seed)) if args.rng_seed else "")
        out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();

        if (args.subst[0] == 'JC'):
            add_arg_col(args.fname, ['alpha'], [args.subst[1]], args.numFamilies)
        elif args.subst[0] == 'KM':          # KM method
            add_arg_col(args.fname, ['alpha', 'beta'], [args.subst[1], args.subst[2]], args.numFamilies)


if __name__=="__main__":
    main()

