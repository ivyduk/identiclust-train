import sys
import re, os
from shutil import copyfile

def tag_fasta(samples, sp, start):

    start = start
    specie = sp
    print (specie)

    IDS = {}

    for i, smp in enumerate(samples):
        ID= specie +"-"+str(i+int(start))
        filename = smp
        tempfile = smp + 'temp'
        o = open(tempfile,"w")
        f = open(filename,'r+')
        for line in f.readlines():
            temp = re.sub(r'>(.+)',r'>'+ID+'|'+'\g<1>', line)
            line = line.replace(line,temp)
            o.write(line)
        o.close()
        IDS[ID]=filename
        copyfile(tempfile, filename)
        os.remove(tempfile)
    return IDS
