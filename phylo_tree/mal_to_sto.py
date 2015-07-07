import re
import sys

"""Change the _whole.mal file to _whole.sto format so that HMMer could use"""

"""The output of readMalFile is an array of lines
ancestor  XXXXX
name1   coord1   coord2   XXXXX
name2   coord1   coord2   XXXXX
"""
def readMalFile(file):
    fileptr = open(file)
    arrForAll = []   
    firstLineFlag = 1  # first line is ancestor
    for line in fileptr:
        line = line.rstrip()
        if firstLineFlag == 1:
            arr = re.split('\s+', line)
            arrForAll.append(arr)
            firstLineFlag = 0
            line = fileptr.next()
        line = line.rstrip()
        arr = re.split('\s+', line)
        new_arr = [arr[0], arr[3], arr[4], arr[5]]
        arrForAll.append(new_arr)
    return arrForAll

"""The input is the output of readMalFile, the output is .afa file"""
def changeToAfa(arr):
    fa = open('new.afa', 'w')
    ances = arr[0][1]   # ancestor seq
    #print ances
    ances_len = len(ances)
    for i in range(1, len(arr)):
        old_line = arr[i]
        new_seq = []
        start_coord = int(arr[i][1])
        num_dash = 0       # count how many dashes needed at the beginning
        count_char = 0     # count character
        for j in range(ances_len):
            if count_char == start_coord + 1:
                #print "num_dash:  ", num_dash
                new_seq = ['~' for x in range(num_dash-1)]
                break
            else:
                num_dash = num_dash + 1
                #print num_dash, ances[j]
                if ances[j] != '*':
                    count_char = count_char + 1
        arr[i][3] = arr[i][3].replace('*', '.').upper().replace('-', '.')
        new_seq.append(arr[i][3])
        new_seq = "".join(new_seq).ljust(ances_len, '~')
        fa.write(">" + arr[i][0][:15] + "\n" + new_seq + "\n")

def changeToSto(arr, output_file):
    fa = open(output_file, 'w')
    fa.write("# STOCKHOLM 1.0" + '\n')
    ances = arr[0][1]   # ancestor seq
    #print ances
    ances_len = len(ances)
    for i in range(1, len(arr)):
        old_line = arr[i]
        new_seq = []
        start_coord = int(arr[i][1])
        num_dash = 0       # count how many dashes needed at the beginning
        count_char = 0     # count character
        for j in range(ances_len):
            if count_char == start_coord + 1:
                #print "num_dash:  ", num_dash
                new_seq = ['.' for x in range(num_dash-1)]
                break
            else:
                num_dash = num_dash + 1
                #print num_dash, ances[j]
                if ances[j] != '*':
                    count_char = count_char + 1
        arr[i][3] = arr[i][3].replace('*', '.').upper().replace('-', '.')
        new_seq.append(arr[i][3])
        new_seq = "".join(new_seq).ljust(ances_len, '.')
        fa.write(arr[i][0][:15] + "." + str(i) + " " + new_seq + "\n")
    fa.write('//')
    fa.close()
    

arr1 = readMalFile(sys.argv[1])
#print arr1
changeToSto(arr1, sys.argv[2])
