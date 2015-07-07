############################################################################
# redhawk.py
# Class for submitting / managing redhawk jobs through a python script
# Written by: John Karro, Jens Muller
# Date: Nov. 1, 2011
#

import re, string, sys
import subprocess
import getopt
import time
import os
import cPickle
from optparse import OptionParser
import pwd

current_user = pwd.getpwuid(os.getuid())[0]  # When run on redhawk with a nohup, os.getlogin() does not work

redhawkStatsRe = re.compile("\s+C\s+[^\s]+\s*$")

class pbsJobHandler:
    """A pbsJobHandler corresponds to a job launched (or to be launched) on redhawk.  Once the object is created (and provided with a command-line execution command),
       the user can extract various inforamtion about the job (current status, output, etc...) and cleanup files."""
    def __init__(self, batch_file, executable, use_pid = True, job_name = None, nodes = 1, ppn = 1, walltime = "40:00:00", address = None, join = False, env = None, queue = None, mail = None, output_location = None, chdir = None, RHmodules = None, file_limit = 6, file_delay = 5): 
        """Constructor.  Requires a file name for the batch file, and the execution command.  Optional parmeters include:
           * use_pid: will embded a process id into the batch file name if true.  Default = true.
           * job_name: A name for the redhawk name.  Default = the batch file name.
           * nodes: number of nodes required for the job.   Default = 1.
           * ppn: number of processors needed for the job.  Default = 1.
           * walltime: Maximum allowed runtime for the job (hours:minutes:seconds).  Default = 40:00:00.  Max. allowed: 400:00:00.
           * mail = when to send email.  Any combination of:
             b   send mail when job begins
             e   send mail when job ends
             a   send mail when job aborts
           * address: additional email addresses to send notification (comma seperated)
           * join: If true, the stdout and stderr files are combined into one file
           * queue: redhawk queue to run on.  Default: redhawk chooses.
           * output_location: Directory to place output files.  Default: current directory.
           * RHmodules: A list of redhawk modules to be loaded before run (e.g. ['Blast+']).  Default: none."""
        self.batch_file_name = batch_file
        if use_pid:
            self.batch_file_name = self.batch_file_name + "." + str(os.getpid())
        self.executable_name = executable
        self.jobname = job_name if job_name else batch_file
        self.join = 'n'
        self.file_limit = file_limit
        self.file_delay = file_delay
        self.status = "unstarted"

        f = open(self.batch_file_name, 'w')

        f.write("#!/bin/bash -l\n")
        s="#PBS -N "+ self.jobname +"\n"
        f.write(s)

     
        #some defaults:
        self.nodes = nodes
        self.ppn = ppn
        self.walltime = walltime
        self.modules = RHmodules
        self.output_location = output_location if output_location else "." 

        s="#PBS -l nodes="+ str(self.nodes)+":ppn="+str(self.ppn)+"\n"
        f.write(s)
        s="#PBS -l walltime="+self.walltime+"\n"
        f.write(s)
        
        if join:
            f.write("#PBS -j oe\n")
            
        if address:    
            s="#PBS -M "+address+"\n"
            f.write(s)

        if queue:
            s="#PBS -q "+queue+"\n"
            f.write(s)

        if env:    
            f.write("#PBS -V\n")

        if mail:
            s="#PBS -m "+mail+"\n"
            f.write(s)

        if output_location:
            s="#PBS -o "+output_location+"\n"
            f.write(s)
            s="#PBS -e "+output_location+"\n"
            f.write(s)

        if chdir:
            s="cd "+chdir+"\n"
            f.write(s)

        else:
            s="cd $PBS_O_WORKDIR\n";
            f.write(s);

        if self.modules != None:
            self.executable_name = "; ".join(["module load " + x for x in self.modules]) + "; " + self.executable_name
        f.write(self.executable_name)

        f.close()
        self.jobid=0;

### submitjob
### Parameters:
###   file is the job script file
###   if preserve is True, don't delete the job script. Delete otherwise
###   optionalFlag is the flag after qsub
###   retry (default set to retry 600 times), the number of times, the job will be submitted in retry
###   seconds between retry (default is 10 seconds)
### return job id if successful
### return -1 if not
    def submitjob(self, preserve=False, print_qsub = False, job_limit = 200, delay=10, user=current_user ):
        """Submit job to redhawk.  Optional parameters:
           * preserve: if False, delete the batch file.  Default = true.
           * job_limit: If the user currently has this many jobs on the batch, wait until one finishes.
           """
        if job_limit > 0:
            self.wait_on_job_limit(limit=job_limit, delay=delay, user=user)

        optionalFlag=''
        retry=600
        sleepTimeBetweenRetry=10
        trial=1;
        cmd = "qsub " + optionalFlag + " " + self.batch_file_name

        while (trial < retry): 
            p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (output, error) = p.communicate()
            if p.returncode == 0:
                break
            trial = trial + 1

        if trial == retry:
            return -1
            
        if not preserve:
            os.remove(self.batch_file_name)
                
        t=re.split('\.',output)
        self.jobid=t[0]
        self.ofile = self.output_location + "/" + self.jobname + ".o" + str(self.jobid)
        self.efile = self.output_location + "/" + self.jobname + ".e" + str(self.jobid)
        if print_qsub:
            print 'qsub jobid', self.jobid
                
        self.status = "running"
        return self

### isJobRunning
### This is primarily useful for waiting a _submitted_ job to finish
###    return False if the job is done, completed for sure
###    return True if the job is in Q, R states [ or that PBS/Torque is not available ]
### Prereq is jobid must be a submitted job
    def isJobRunning ( self ):
        """Query of the object represented by the job is still running."""
        #cmd = "qstat " + str(self.jobid)
        
        #magicString='Unknown Job Id'  ### magicString _might_ need to be changed if Torque version changes
        #(output, error) = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

        if self.ofile_exists():  #output.find(magicString) >=0 or redhawkStatsRe.search(output):
            self.status = "finished"
            return False
        

        return True
    
    def wait_on_job(self, delay=10):
        """Spin until job completes."""
        while self.isJobRunning() == True:
            time.sleep(delay)
        return self.ofile_exists()    

    def ofile_name(self):
        """Get the name of the file containing the job stdio output."""
        return self.ofile

    def efile_name(self):
        """Get the name of the file containing the job stderr output."""
        return self.efile

    def ofile_exists(self):
        """Does the file contiining the job stdout output exist?"""
        return os.path.isfile(self.ofile)

    def efile_exists(self):
        """Does the file contiining the job stderr output exist?"""
        return os.path.isfile(self.efile)


    def ofile_handle(self):
        """Return a handle to the file containing the job stdout output."""
        if not self.status == "finished":
            raise NameError("redhawk: unfinished ofile check")
        tries = 0
        while not self.ofile_exists() and tries < self.file_limit:
            time.sleep(self.file_delay)
            tries = tries+1
        
        if os.path.isfile(self.ofile_name()):
            return open(self.ofile_name(), "r")

        raise NameError("redhawk: unfound ofile")

    def efile_handle(self):
        """Return a handle to the file containing the job stderr output."""
        if not self.status == "finished":
            raise NameError("redhawk: unfinished efile check")

        tries = 0
        while not self.efile_exists() and ties < self.file_limit:
            time.sleep(self.file_delay)
            tries = tries+1
        
        if os.path.isfile(self.efile_name()):
            return open(self.efile_name(), "r")

        raise NameError("redhawk: unfinished efile check")


    def ofile_string(self):
        """Return the entire contents of the stdout file as a single string."""
        fp = self.ofile_handle()
        if (fp):
            return "\n".join([line.rstrip() for line in fp])
        return None

    def efile_string(self):
        """Return the entire contents of the stderr file as a single string."""
        fp = self.efile_handle()
        if (fp):
            return "\n".join([line.rstrip() for line in fp])
        return None

    def erase_files(self):
        """Erase the stdio and stderr files."""
        self.ofile_handle()
        self.efile_handle()

        os.remove(self.ofile_name())
        os.remove(self.efile_name())
        return None
    
    def getResults(self, cleanup=True):
        """Retrieve strings and cleanup"""
        self.wait_on_job()
        stdout_str = self.ofile_string()
        stderr_str = self.efile_string()
        if cleanup:
            self.erase_files()
        return (stdout_str, stderr_str)

    def wait_on_job_limit(self,limit=200,delay=10,user=current_user):
        while 1==1:
            numJobsInQueue = get_number_of_jobs_in_queue(current_user)

            if (numJobsInQueue < limit):
                return None
            time.sleep(delay)


### return the number of jobs in queue, whatever the state is
def get_number_of_jobs_in_queue(user=current_user):
    cmd = "qstat -u "+user + " 2>/dev/null | grep " + user

    output,error = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    return len([True for line in output.split("\n") if line and not redhawkStatsRe.search(line)])
                                                        
         

### The following functions are not well tested -- use with CAUTION.
def storePBS(pbsList, fp):
    for p in pbsList:
        cPickle.dump(pbsList, fp)

def loadPBS(fp):
    return [obj for obj in cPickle.load(fp)]



################################
### Sample code
### exe = "ls *"     # The command we want to run (to run multiple commands, seperate with semi-colons)
### o = pbsJobHandler(batch_file = "batch.txt", executable = exe, mail = "bea");    # Set up the job
### o.submitjob()                        # Submit the job to a redhawk queue
### if o.isJobRunning(): print "yes"     # Check to see if job is still in the queue or running
### o.wait_on_job()                      # "spin" until job is finished
### output_string = o.ofile_string()     # Get the ouput
### o.erase_files()                      # Erase the output files
