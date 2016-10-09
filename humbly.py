#!/usr/bin/python
#------------------------------------------------------------------------------
# Filename: humbly.py
# Created: 06.10.2016
# Author: Milan Kovacevic, kovmil@gmail.com
#------------------------------------------------------------------------------
import sys
import re
from collections import Counter
from collections import namedtuple

#------------------------------- Functions ------------------------------------
#Function for parsing string (line from .pileup file) into PileupStruct structure
def str_to_pileup_struct(str):
    spl = str.split()
    #Checking if coverage is >0, otherwise error occur
    if int(spl[3]) > 0:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=spl[4], quality=spl[5])
    else:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=0, quality=0)
    return pile;

#Function for determining quality of chosen ALT based on quality string provided
def quality(base, qual, alt):
    positions = [m.start() for m in re.finditer(str(alt), base.replace(",",".").upper())]
    num_of_chosen = 1
    if len(positions) > 0:
        num_of_chosen = len(positions)
        s = 0
        for n in positions:
            s += ord(qual[n]) #ord(c) returns ASCII value of 'c'
        return round(s*1.0/num_of_chosen, 2);

#------------------------------- Main program ---------------------------------
if len(sys.argv) == 1:
    print "Error: Usage: ./humbly.py file_name.pileup"
    exit()

if sys.argv[1][-6:] != "pileup":
    print "Error: Please provide PILEUP file format"
    exit()

with open(sys.argv[1]) as f:
    lines = f.readlines()

#Variable for printing VCF file
vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tGENOTYPE"

#C-like Typedef struct PiledupStruct
PileupStruct = namedtuple("PileupStruct", "chrom position ref coverage bases quality")

prog = re.compile(r'(\.*[ATCGN]+\.*)+')

for iter in lines:
    piled_up = str_to_pileup_struct(iter)

    alt = None
    genotype = "0/0"

    #TODO:Quality of the base string and ref. genome will be added in later versions
    #Regex for finding bases which are considered
    #Ignoring DNA direction (every ',' is '.', 'a' is 'A' etc.)
    result = prog.match(str(piled_up.bases).replace(",",".").upper())
    if result != None:

        #Counter is dict subclass for counting the number of occurances of every character in string
        x = Counter(result.group())
        mc = x.most_common(1)[0][0]

        #Case when there is same number of different bases (e.g. AAAACCCC)
        newdict = {}
        for k,v in x.iteritems():
            newdict.setdefault(v, []).append(k)
        for k,v in newdict.iteritems():
            if len(v)==2:
                if v[1]!='.':
                    alt = v #Genotype 1/2
                    genotype = "1/2"
                else:
                    alt = v[0] #Genotype 0/1
                    genotype = "0/1"

        #Case when second_most_common base is still considerable
        if len(x.keys()) > 1:
            smc = sorted(x.items(), key=lambda x:x[1])[0][0]
            smc_cnt = sorted(x.items(), key=lambda x:x[1])[0][1]
        if (alt == None and smc != '.' and len(result.group())*1.0/4 <= smc_cnt and len(x.keys()) > 1):
            alt = smc
            print str(result.group())
            genotype = "0/1"

        if (mc != '.' and alt == None):
            alt = mc; #Genotype 1/1
            genotype = "1/1"

    q = quality(str(piled_up.bases),str(piled_up.quality), alt)

    if alt != None:
        vcf += "\n"+str(piled_up.chrom)
        vcf += "\t"+str(piled_up.position)
        vcf += "\t"+str(".")
        vcf += "\t"+str(piled_up.ref)
        vcf += "\t"+str(alt)
        vcf += "\t"+str(q)
        vcf += "\t"+str(genotype)

print vcf
