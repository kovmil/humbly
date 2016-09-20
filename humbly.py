#!/usr/bin/python
import sys
import re
from collections import Counter
from collections import namedtuple

#------------------------------- Functions ---------------------------------
#Function for parsing string (lane from .pileup file) into PileupStruct structure
def str_to_pileup_struct(str):
    spl = str.split()
    #Checking if coverage is >0, otherwise error occur
    if int(spl[3]) > 0:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=spl[4], quality=spl[5])
    else:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=0, quality=0)
    return pile;

#Function for picking ALT base using string of bases and ref. genome
#Quality of the base string will be added in later versions
def alt_pick(base):
    #Regex for finding bases which are consi
    prog = re.compile(r'(\.*[ATCGNatcgn]+\.*)+')
    result = prog.match(base)
    if result is None:
        return None

    #Counter is dict subclass for counting the number of occurances of every character in string
    x = Counter(result.group())

    #Case when there is same number of different bases (ie. AAAACCCC)
    newdict = {}
    for k,v in x.iteritems():
        newdict.setdefault(v, []).append(k)
    for k,v in newdict.iteritems():
        if len(v)==2:
            if v[1]!='.':
                return v #Genotype 1/2
            else:
                return v[0] #Genotype 1/1

    if (x.most_common(1)[0][0] != '.'):
        return x.most_common(1)[0][0]; #Genotype 0/1

#Function for determining quality of chosen ALT based on quality string provided
def quality(base, qual, alt):
    positions = [m.start() for m in re.finditer(str(alt), base.replace(",",".").upper())]
    num_of_chosen = 1
    if len(positions) > 0:
        num_of_chosen = len(positions)
        s = 0
        for n in positions:
            s += ord(qual[n]) #ord(c) returns ASCII value of 'c'
        return round(s*1.0/num_of_chosen, 5);

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
vcf = []
vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL"

#C-like Typedef struct PiledupStruct
PileupStruct = namedtuple("PileupStruct", "chrom position ref coverage bases quality")

for iter in lines:
    piled_up = str_to_pileup_struct(iter)

    #Ignoring DNA direction (every ',' is '.', 'a' is 'A' etc.)
    alt = alt_pick(str(piled_up.bases).replace(",",".").upper())
    
    q = quality(str(piled_up.bases),str(piled_up.quality), alt)

    if alt != None:
        vcf += "\n"+str(piled_up.chrom)
        vcf += "\t"+str(piled_up.position)
        vcf += "\t"+str(".")
        vcf += "\t"+str(piled_up.ref)
        vcf += "\t"+str(alt)
        vcf += "\t"+str(q)

print vcf
