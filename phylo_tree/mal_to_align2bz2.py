import re
import sys
import bz2
import os

"""This program converts the .mal file to .align2 file so that TEprofiler.py can use it"""

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

"""The input is the output of readMalFile, the output is .align2 file"""
def changeToAlign2(arr):
    #fa = open('sim.align2', 'w')
    fa = bz2.BZ2File(proc_id + "_sim.align2.bz2", 'w')
    ances = arr[0][1]   # ancestor seq
    for i in range(1, len(arr)):
        old_line = arr[i]
        fa.write(arr[i][0] + " " + proc_id + "_SimuRepBase" + " " + arr[i][1] + " " + arr[i][2] + "\n")
        organism = []
        original = []
        seq = arr[i][3]
        insert_counter = 0
        for j in range(len(seq)):
            if seq[j].isupper() and seq[j]!='-':      # match
                organism.append(seq[j])
                original.append(ances[int(arr[i][1])+j-insert_counter])
            elif seq[j].isupper()==False and seq[j]!='-':   # insert
                organism.append(seq[j].upper())
                original.append('-')
                insert_counter = insert_counter + 1
            elif seq[j]=='-':                         # delete
                organism.append('-')
                original.append(ances[int(arr[i][1])+j-insert_counter])
        fa.write("".join(organism) + '\n' + "".join(original) + '\n')
    fa.close()


proc_id = sys.argv[1]
arr = readMalFile(proc_id + "_sim_out.mal")
changeToAlign2(arr)
                
                
