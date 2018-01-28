import sys
import re, os
from shutil import copyfile
import sys, getopt, json

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
    print (IDS)

    with open("dics/%s.json" %sp, 'w') as f:
        json_str = json.dumps(IDS, indent=4)
        f.write(json_str)
    f.close()
    return IDS
