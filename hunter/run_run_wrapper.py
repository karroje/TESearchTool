import sys
import subprocess
import re
import os
import fileinput
import argparse

def main():
    
    cmd = "python run_wrapper.py -ph -fn result.txt"
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();

    for i in range(20, 31):
        #print "run_run_wrappyer: i=", i
        #cmd = "python run_wrapper.py -s 0.03 -fn testResult.txt --seed " + str(i)
        #cmd = "python run_wrapper.py -s KM 0.01 0.02 -fn testResult.txt"
        try:    
           cmd = "python run_wrapper.py -s JC " + str(i*0.01) +" -nf 5 -fn result.txt -rd EXPO 5.0 5.0 --seed " + str(i)
           out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();
        except:
           print "Unexpected error: ", sys.exc_info()[0]

if __name__=="__main__":
    main()
