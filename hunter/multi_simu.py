import sys
import subprocess
import os
import re

def simulator():
    cmd = "python basic_simulator_multi.py sim_params.txt sim_out.mal sim_out.fa true_coords.txt; python mal_to_align2bz2.py"
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();

def pipeline():
    cmd = "python pipeline.py SimuRepBase SimuRepBase.fa sim_out.fa 12 hmmsearch_sim.output sim.align2.bz2"
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();

def analyser():
    cmd = "python result_analyser.py"
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();
    arr = re.split("\n", out)
    print arr

if __name__=="__main__":
    simulator()
    pipeline()
    analyser()
