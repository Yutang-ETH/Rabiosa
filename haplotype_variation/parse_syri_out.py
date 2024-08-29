#!/usr/bin/env python3

# parse the syri.out by python
# want to output three files, SNP, small INDEL (size < 50), large INDEL (size >= 50)

import sys

myinput = sys.argv[1]
SNPoutput = sys.argv[2]
INDELoutput = sys.argv[3]
PAVoutput = sys.argv[4]

mySNP = []
myINDEL = []
myPAV = []

with open(myinput, "r") as f:
    for line in f:
        if not line.startswith("#"):
            if line.split()[10] == "SNP":
                chr = line.split()[0]
                start = line.split()[1]
                end = line.split()[2]
                length = "1"
                type = line.split()[10]
                ref = line.split()[3]
                alt = line.split()[4]
                myline = "\t".join([chr, start, end, ref, alt, str(length), type])
                mySNP.append(myline)
            elif line.split()[10] == "INS":
                if len(line.split()[4]) - len(line.split()[3]) < 50:
                    chr = line.split()[0]
                    start = line.split()[1]
                    end = line.split()[2]
                    length = len(line.split()[4]) - len(line.split()[3])
                    ref = line.split()[3]
                    alt = line.split()[4]
                    type = "INDEL"
                    myline = "\t".join([chr, start, end, ref, alt, str(length), type])
                    myINDEL.append(myline)
                else:
                    chr = line.split()[0]
                    start = line.split()[1]
                    end = line.split()[2]
                    length = len(line.split()[4]) - len(line.split()[3])
                    ref = line.split()[3]
                    alt = line.split()[4]
                    type = "INS"
                    myline = "\t".join([chr, start, end, ref, alt, str(length), type])
                    myPAV.append(myline)
            elif line.split()[10] == "DEL":
                if len(line.split()[3]) - len(line.split()[4]) < 50:
                    chr = line.split()[0]
                    start = line.split()[1]
                    end = line.split()[2]
                    length = len(line.split()[3]) - len(line.split()[4])
                    ref = line.split()[3]
                    alt = line.split()[4]
                    type = "INDEL"
                    myline = "\t".join([chr, start, end, ref, alt, str(length), type])
                    myINDEL.append(myline)
                else:
                    chr = line.split()[0]
                    start = line.split()[1]
                    end = line.split()[2]
                    length = len(line.split()[3]) - len(line.split()[4])
                    ref = line.split()[3]
                    alt = line.split()[4]
                    type = "DEL"
                    myline = "\t".join([chr, start, end, ref, alt, str(length), type])
                    myPAV.append(myline)

with open(SNPoutput, "w") as SNPf:
    for i in range(0, len(mySNP)):
        SNPf.write(mySNP[i] + "\n")

with open(INDELoutput, "w") as INDELf:
    for i in range(0, len(myINDEL)):
        INDELf.write(myINDEL[i] + "\n")

with open(PAVoutput, "w") as PAVf:
    for i in range(0, len(myPAV)):
        PAVf.write(myPAV[i] + "\n")




             
