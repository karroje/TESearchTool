import re
import sys
from Bio.Seq import Seq


# KARRO: Note that I replaced "human" with "organism" in all places.  This isn't
# going to be used only on humans, so calling it human makes the code confusing.


def findAncestorSeq(file, family):
    fileptr1 = open(file)
    
    correctFamily1 = False;          #this part is to get the family from .fa and make it one whole line   
    temparr = []                     #to store the whole line
    r = re.compile(family)
    for line in fileptr1:
        line = line.rstrip()

        if line[0] == ">":
            if correctFamily1:
                break;
            if (r.search(line)):
                correctFamily1 = True
                temparr.append(line)
            else:
                correctFamily1 = False;
        else:
            if correctFamily1:
                temparr.append(line)
                #print line
    fileptr1.close()

    #temparr = temparr[0] + "\n" + "".join(temparr[1:])
    return "".join(temparr[1:])
    #print temparr;    



def readFamilySequences(file, temparr, family):
    fileptr2 = open(file);

    correctFamily2 = 0;      #this part is to find the family in .align2 and reduce dashes
    arrforall = []    #the whole arr
    arrelement = []   #each element of arrforall contains (origseq,organismseq,coordinates in MIRb)
    coor_arr = []     #the coordinates in MIRb line
    for line in fileptr2:
        if line[0] == '#':    # KARRO: Allow us to comment out lines int the file (for testing)
            continue
        line = line.rstrip()
        arr = re.split("\s+", line)
        if len(arr) > 1:
            if arr[1] == family:
                coor_arr.append(int(arr[2]))  # KARRO: Changed these to ints for consistancy. 
                coor_arr.append(int(arr[3]))
                correctFamily2 = 1
                line = "\t".join(arr)
            else:
                correctFamily2 = 0;
        else:
            if correctFamily2 == 1:        #organism sequece line
                correctFamily2 = 2
                organism_sequence = line
            elif correctFamily2 == 2:      #original sequnce line
                origi_sequence = line
                line = line.lower()
                arr2 = re.split("-+", line)    #reduce dashes in "line"
                line = "".join(arr2)
                r0 = re.search(line, temparr)  #search modified origi in MIRb
                if r0 is None:
                    line = Seq(line)       #to use biopython
                    line = line.reverse_complement()
                    line = str(line)
                    organism_sequence = Seq(organism_sequence)
                    organism_sequence = organism_sequence.reverse_complement()
                    organism_sequence = str(organism_sequence)
                    origi_sequence = Seq(origi_sequence)
                    origi_sequence = origi_sequence.reverse_complement()
                    origi_sequence = str(origi_sequence)
                    
                temparr2 = temparr.replace(line, origi_sequence)  #change MIRb to original
                #print temparr2        #MIRb with original changed
                #print line            #lower without dashes original piece
                pat = origi_sequence
                r1 = re.search(pat, temparr2)
                start = r1.start()
                end = r1.end()
                arrelement = [origi_sequence, organism_sequence, start, start, coor_arr[0], coor_arr[1]]
                arrforall.append(arrelement)
                coor_arr = []
                #print arrforall
                #print "\n"
    return arrforall

def extract(family, organismFile, ancestorFile, print_aligned, insert_char, chr_num = "", print_insert = False):
    fileptr2 = open(organismFile)
    if not fileptr2:
        print "Bad file";
        return;

    # KARRO: Replaced "temparr" with "ancestorSeq" -- a more correct description

    # Get the ancestor sequence
    ancestorSeq = findAncestorSeq(ancestorFile, family)

    # Get the information for the modern sequence, return as a list of lists.  Each sublists contains:
    # 1) Ancestor portion of the alignment string
    # 2) Decendant portion of the alignment string
    # 3) Start position of the alignment in the ancestor sequence relative to the unmodified string
    # 4) Start position of the alginment in the ancestor sequence relative to the *modified* ancestor string (where dashes have been aded)
    # 5) Start position of the modern TE relative to its chromosome
    # 6) End position of the modern TE relative to its chromosome
    alignmentRecords = readFamilySequences(organismFile, ancestorSeq, family);


    for j in range(len(alignmentRecords)):
        (ori, organism, start, realstart, coor1, coor2) = alignmentRecords[j]
        counter1 = realstart
        counter2 = 0           #coor of ori and hum


        while (counter2 < len(organism)):
            if ancestorSeq[counter1] != insert_char and ori[counter2] != insert_char and ori[counter2] != "-":
                counter1 = counter1 + 1
                counter2 = counter2 + 1
            elif ancestorSeq[counter1] != insert_char and ori[counter2] == "-":   #need to add a dash in MIRb
                ancestorSeq = ancestorSeq[:counter1] + insert_char + ancestorSeq[counter1:]
                for i in range(len(alignmentRecords)):   #modify all the genes except the one in the outer forloop
                    (ori2, organism2, start2, realstart2, coor_1, coor_2) = alignmentRecords[i]
                    if organism2 != organism:   
                        if realstart2 < counter1 and realstart2 + len(organism2) > counter1 and i<j: #add a dash in ori and hum if in that scope
                            ori2 = ori2[:(counter1 - realstart2)] + insert_char + ori2[(counter1 - realstart2):]
                            organism2 = organism2[:(counter1 - realstart2)] + insert_char + organism2[(counter1 - realstart2):]
                        if realstart2 >= counter1:   #in this case only change the coordinates of the genes
                            realstart2 = realstart2 + 1
                    alignmentRecords[i] = (ori2, organism2, start2, realstart2, coor_1, coor_2)
                counter1 = counter1 + 1
                counter2 = counter2 + 1
            elif ancestorSeq[counter1] == insert_char and ori[counter2] != insert_char and ori[counter2] != "-":   #need to add a dash in ori and hum
                if counter2 == 0:
                    counter1 = counter1 + 1
                    realstart = realstart + 1
                    alignmentRecords[j] = (ori, organism, start, realstart, coor1, coor2)
                else:
                    ori = ori[:counter2] + insert_char + ori[counter2:]
                    organism = organism[:counter2] + insert_char + organism[counter2:]    
                    alignmentRecords[j] = (ori, organism, start, realstart, coor1, coor2)
                    counter1 = counter1 + 1
                    counter2 = counter2 + 1
            #elif ancestorSeq[counter1] == "-" and ori[counter2] == "-":
            elif (ancestorSeq[counter1] == "-" or ancestorSeq[counter1] == insert_char) and (ori[counter2] == "-" or ori[counter2] == insert_char):
                counter1 = counter1 + 1
                counter2 = counter2 + 1
    
    for i in range(len(alignmentRecords)):                      #make organism lower case corresponding to dashes in ancestor
        (ori, organism, start, realstart, coor1, coor2) = alignmentRecords[i]
        organismlist = []
        organismlist = list(organism)
        for j in range(len(ori)):
            if ori[j] == "-" and organism[j] != "-":
                organismlist[j] = organismlist[j].lower()
        organismlist = "".join(organismlist)
        alignmentRecords[i] = (ori, organismlist, start, realstart, coor1, coor2)  #use organismlist (in fact it's string now) to replace organism


    #extract chromosome name
    if not chr_num:
        chr_num = organismFile[0:organismFile.find('.')]

    if print_aligned:
        sys.stdout.write("{0:10}\t{1:10}\t{2:10}\t{3:10}\t{4:10}\t{5}".format("ancestor", "", "", "", "", ancestorSeq) + "\n")
    else:
        if not print_insert:
            ancestorSeq = ancestorSeq.replace(insert_char, "");
        sys.stdout.write("ancestor\t" + ancestorSeq + "\n")

    for (ori, organism, start, realstart, coor1, coor2) in alignmentRecords:    #print out all the organism genes and coordinates
        end = start
        for j in range(0, len(ori)):
            if ori[j] != "-" and ori[j] != insert_char:
                end = end + 1

        if print_aligned:
            sys.stdout.write("{0:10}\t{1:10}\t{2:10}\t{3:10}\t{4:10}\t".format(chr_num, coor1, coor2, start, end))
            l = len(re.match("(\%s*[^\%s]){%d}" % (insert_char, insert_char, start+1), ancestorSeq).group(0)) - 1
            sys.stdout.write(" "*l + organism + "\n")
        else:
            if not print_insert:
                organism = organism.replace(insert_char, "");
            sys.stdout.write("\t".join([str(x) for x in [chr_num, coor1, coor2, start, end, organism]]) + "\n")
    
if __name__ == "__main__":
    if len(sys.argv)==1:
        extract("MIRb", "chrT.fa.align2", "test.repbase.fa", False, '*', "")
    else:
        marker = 1;
        print_aligned = False;
        print_insert = False
        insert_char = '*'
        chr_num = ""    # If empty, will extract the chromsome number from the organism file name, assumed form "chrX."

        while sys.argv[marker][0] == '-':
            switch = sys.argv[marker]
            if switch == '-a':
                print_aligned = True
            elif switch == '-i':
                marker = marker + 1
                insert_char = sys.argv[marker]
            elif switch == '-c':
                marker = marker + 1
                chr_num = sys.argv[marker]
            elif switch == '-pi':
                print_insert = True;
            marker = marker + 1

        family, chr_file, ancestor_file = sys.argv[marker:]

        extract(family, chr_file, ancestor_file, print_aligned, insert_char, chr_num, print_insert);

else:
    #extract("MIRb", "chrT.fa.align2", "test.repbase.fa", True, insert_char)
    extract("MIRb", "2chr22.fa.align2", "2Homo_sapiens.repbase.fa", "version1")
    

