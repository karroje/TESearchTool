#!/usr/bin/python

import sys
import subprocess




def main(numtrials):


    
    #proc = subprocess.Popen(['mkdir', '.temp'], shell=False)
  
    # for each data entry 
    #   generate data ./simulation/basic_simulator.pl <sim_params.txt>
    #   build hash ./translate <input> <output> <anch_len>
    #   wait a year
    #   run tester tester -p <hmm> -i <info> -h <hash>
    #   output results
    #   delete .hash, .fa, .hmm, .hmm.info, .mal
    
    ANCHOR_SIZE = 12
    for i in range(numtrials):
        print ("Trial#" + str(i))

        p1 = subprocess.Popen(['./basic_simmulator.pl', 
                              'sim_params.txt',
                              './.temp/trial' + str(i) + '.mal',
                              './.temp/trial' + str(i) + '.fa'],
                              shell=False,
                              stdout=subprocess.PIPE)

        p1.wait()

        p3 = subprocess.Popen(['../build_hmm',
                              './.temp/trial' + str(i),
                              './.temp/trial' + str(i) + '.mal'],
                              shell=False,stdout=subprocess.PIPE)
        p3.wait()

        p4 = subprocess.Popen(['../tester',
                               '-p',
                               './.temp/trial' + str(i) + '.hmm',
                               '-i',
                               './.temp/trial' + str(i) + '.hmm.info',
                               '-h',
                               './.temp/trial' + str(i) + '.hash',
                               '-q',
                               './.temp/trial' + str(i) + '.fa'],
                              shell=False,
                              stdout = subprocess.PIPE)
        p4.wait()
        out = p4.communicate()[0]
        print(out)
                         

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Bad Params"
        sys.exit()

    num_trials = int(sys.argv[1])
    main(num_trials)
