import sys
import subprocess
import re
import os

#proc_id = str(os.getpid())

def simulator(seed):
    cmd_str_part2 = ""
    for i in range(int(num_family)):
        cmd_str_part2 = cmd_str_part2 + "; python mal_to_align2bz2.py " + str(i) + "_" + proc_id
    cmd = "python basic_simulator_multi.py " + (("--seed " + seed) if seed else "") + " sim_params.txt sim_out.mal " \
          + "sim_out.fa true_coords.txt " + proc_id + " " + num_family + cmd_str_part2
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();

def pipeline():
    cmd = "python pipeline.py " + proc_id + "_SimuRepBase " + proc_id + "_SimuRepBase.fa " \
          + proc_id + "_sim_out.fa 12 " + proc_id + "_hmmsearch_sim.output " + proc_id + " " + num_family + " " + \
		  proc_id + "_sim.align2.bz2"
    out, err = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();

def analyser():
    cmd = "python result_analyser.py %d_" + proc_id
    output_arr_list = []                              # each "arr" is a result output of one simulated family
    for i in range(int(num_family)):
        out, err = subprocess.Popen(cmd%(i), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate();
        arr = re.split("\n", out)
        output_arr_list.append(arr)

    #print "Opening: ", result_file_name
    fa = open(result_file_name, 'a')
    for arr in output_arr_list:                       # for the result of each family
        for i in range(len(arr)):                     # for each line in the result
            sub_arr = re.split("\s+", arr[i])
            #print str(sub_arr)
            if i==0:
                fa.write(sub_arr[2].ljust(10)+sub_arr[5].ljust(10)+sub_arr[8].ljust(10)+sub_arr[11].ljust(10))
            else:
                fa.write(sub_arr[-1].ljust(20))
        fa.write('\n')
    fa.close()
    

def main(seed):
    simulator(seed)
    pipeline()
    analyser()


if __name__=="__main__":
    proc_id = sys.argv[1]
    result_file_name = sys.argv[2]
    num_family = sys.argv[3]
    seed = sys.argv[4] if len(sys.argv) > 4 else None;
    print "one_whole_sim.py: ", result_file_name, seed
    main(seed);
