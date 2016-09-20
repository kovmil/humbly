#!/usr/bin/python
import sys
import re
from collections import Counter
from collections import namedtuple
#typedef struct PiledupStruct
PileupStruct = namedtuple("PileupStruct", "chrom position ref coverage bases quality")

#function for parsing string into PileupStruct structure
def str_to_pileup_struct(str):
    spl = str.split()
    #checking if coverage is >0
    if int(spl[3]) > 0:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=spl[4], quality=spl[5])
    else:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=0, quality=0)
    return pile;

#function for picking ALT base (quality not included in calculations)
def alt_pick(base):
    #regex for finding ..A.A...A..A.AAAC.C.,,. strings
    prog = re.compile(r'(\.*[ATCGNatcgn]+\.*)+')
    result = prog.match(base)
    if result is None:
        return None

    #Counter is counting nuber of occurances of every character in string
    x = Counter(result.group())
    #most common character is chosen (quality will be added in calculations)

    #case when there is same number of different bases (ie. AAAACCCCG..)
    newdict = {}
    for k,v in x.iteritems():
        newdict.setdefault(v, []).append(k)
    for k,v in newdict.iteritems():
        if len(v)==2:
            #print v
            if v[1]!='.':
                return v #1/2
            else:
                return v[0] #1/1

    if (x.most_common(1)[0][0] != '.'):
        return x.most_common(1)[0][0]; #0/1


#function for determining quality based on quality string passed of bases and ALT chosen
def quality(base, qual, alt):
    positions = [m.start() for m in re.finditer(str(alt), base.replace(",",".").upper())]
    num_of_chosen = 1
    if len(positions) > 0:
        num_of_chosen = len(positions)
        s = 0
        for n in positions:
            s += ord(qual[n]) #ord(c) returns ASCII value of 'c'
        return round(s*1.0/num_of_chosen, 5);

#------------------------------- main program ---------------------------------
#number of arguments checking
if len(sys.argv) == 1:
    print "Error: Usage: ./humbly.py file_name.pileup"
    exit()

if sys.argv[1][-6:] != "pileup":
    print "Error: Please provide PILEUP file format"
    exit()

#opening file and iterating line by line
with open(sys.argv[1]) as f:
    lines = f.readlines()

#vcf output - list of strings
vcf = []
vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL"

for iter in lines:
    piled_up = str_to_pileup_struct(iter)
    #ignoring DNA direction (every ',' is '.' and 'a' is 'A' etc.)
    alt = alt_pick(str(piled_up.bases).replace(",",".").upper())
    q = quality(str(piled_up.bases),str(piled_up.quality), alt)

    if alt != None:
        vcf += "\n"+str(piled_up.chrom)
        vcf += "\t"+str(piled_up.position)
        vcf += "\t"+str(".")
        vcf += "\t"+str(piled_up.ref)
        vcf += "\t"+str(alt)
        vcf += "\t"+str(q)
        #print "Quality:"+str(q)+" Bases:"+str(piled_up.bases.replace(",",".").upper())+" ALT:"+str(alt)
print vcf
