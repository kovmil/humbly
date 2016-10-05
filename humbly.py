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
    prog = re.compile(r'(\.*,*[ATCGNatcgn]+\.*,*)+')
    result = prog.match(base)
    if result is None:
        return None
    #ignoring DNA direction (every ',' is '.' and 'a' is 'A')
    #Counter is counting nuber of occurances of every character in string
    x = Counter(result.group().replace(",",".").upper())
    #most common character is chosen (quality will be added in calculations)
    return x.most_common(1)[0][0];

#function for determining quality based on quality string passed, bases and ALT chosen
def quality(base, qual, alt):
    positions = [m.start() for m in re.finditer(str(alt), base.replace(",",".").upper())]
    return positions;

#main program
with open('example.pileup') as f:
    lines = f.readlines()
#vcf output - list of strings
vcf = []
for iter in lines:
    piled_up = str_to_pileup_struct(iter)
    alt = alt_pick(str(piled_up.bases))
    q = quality(str(piled_up.bases),str(piled_up.quality), alt)

    if alt != None and alt != '.':
        print "Positions: "+str(q)+" Bases: "+str(piled_up.bases.replace(",",".").upper())+" Chosen base: "+str(alt)
    #    print alt
