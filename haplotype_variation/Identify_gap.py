#!/usr/bin/env python3

# I want to find gaps in an assembly and output a bed file to show the position of the gap

import sys
import re

myinput = sys.argv[1]
myoutput = sys.argv[2]

# import the fasta file

myheader = []
mysequence = []

with open(myinput, "r") as f:
    for line in f:
        if line.startswith(">"):
            myheader.append(line.strip(">").strip())
        else:
            mysequence.append(line.strip())

# check the length of myheader and mysequence
print(len(myheader), len(mysequence))

# find gaps in sequence

with open(myoutput, "w") as outf:
    for i in range(len(myheader)):
        for x in re.finditer("N+", mysequence[i]):
            outf.write(("\t").join([myheader[i], str(x.start()), str(x.end()), "\n"]))

