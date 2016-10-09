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

SIGNIFICANT_PERCENTAGE_OF_BASES = (0.25)
MAJOR_PERCENTAGE_OF_BASES = (0.75)

vcf = ""

#------------------------------- Functions ------------------------------------
#Function for parsing string (line from .pileup file) into PileupStruct structure
def str_to_pileup_struct(str):
    spl = str.split()
    #Checking if coverage is >0, otherwise error occur
    if int(spl[3]) > 0:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=spl[4].replace(",",".").upper(), quality=spl[5])
    else:
        pile = PileupStruct(chrom=spl[0], position=spl[1], ref=spl[2], coverage=spl[3], bases=0, quality=0)
    return pile;

#Function for determining quality of chosen ALT based on quality string provided
def quality(base, qual, alt):
    positions = [m.start() for m in re.finditer(str(alt), base)]
    num_of_chosen = 1
    if len(positions) > 0:
        num_of_chosen = len(positions)
        s = 0
        for n in positions:
            s += ord(qual[n]) #ord(c) returns ASCII value of 'c'
        return round(s*1.0/num_of_chosen, 2);

def header():
    global vcf
    vcf += "##fileformat=VCFv4.2\n"
    vcf += "##FILTER=<ID=PASS, Description=\"All filters passed\">\n"
    vcf += "##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n"
    vcf += "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n"
    vcf += "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"
    vcf += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"



    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tGENOTYPE\t" + str(sys.argv[1]).replace(".pileup", ".bam")

#------------------------------- Main program ---------------------------------
if len(sys.argv) == 1:
    print "Error: Usage: ./humbly.py file_name.pileup"
    exit()

if sys.argv[1][-6:] != "pileup":
    print "Error: Please provide PILEUP file format"
    exit()

with open(sys.argv[1]) as pileup_file:
    lines = pileup_file.readlines()

#Variable for printing VCF file
header()


#C-like Typedef struct PiledupStruct
PileupStruct = namedtuple("PileupStruct", "chrom position ref coverage bases quality")

reg_exp = re.compile(r'(\.*((\+|-)[0-9]+[ATCGN]+)*[ATCGN.]+\.*)+')
indel_exp = re.compile(r'(\.*[\^$]*((\+|-)[0-9]+[ATCGN]+)*[ATCGN.]+\.*)+')



for iter in lines:
    piled_up = str_to_pileup_struct(iter)

    alt = None
    genotype = "0/0"

    #Regex for finding bases which are considered
    #Ignoring DNA direction (every ',' is '.', 'a' is 'A' etc.)
    match_found = reg_exp.match(str(piled_up.bases))
    indel_match = indel_exp.match(str(piled_up.bases))
    if match_found != None:

        #Counter is dict subclass for counting the number of occurances of every character in string
        base_count = Counter(str(match_found.group()))
        #if base_count['+'] > 0 or base_count['-'] > 0:
            #print indel_match.group()

        mc_list = base_count.most_common(2)
        base_len = len(str(match_found.group()))

        mc = mc_list[0][0]
        mc_cnt = mc_list[0][1]

        if len(mc_list) > 1:
            smc = mc_list[1][0]
            smc_cnt = mc_list[1][1]
        else:
            smc = None
            smc_cnt = 0

        if mc == '.':
            if smc_cnt > base_len/4:
                alt = smc
                genotype = "0/1"
            else:
                alt = piled_up.ref
                genotype = "0/0"
        else:
            alt = mc
            if mc_cnt > 3*base_len/4:
                genotype = "1/1"
            elif smc_cnt > base_len/4:
                if smc == '.':
                    genotype = "0/1"
                else:
                    alt = [alt, smc]
                    genotype = "1/2"

    if genotype != "0/0":
        q = quality(str(piled_up.bases),str(piled_up.quality), alt)

        vcf += "\n"+str(piled_up.chrom)
        vcf += "\t"+str(piled_up.position)
        vcf += "\t"+str(".")
        vcf += "\t"+str(piled_up.ref)
        vcf += "\t"+str(alt)
        vcf += "\t"+str(q)
        vcf += "\t"+str(genotype)

print vcf
#TODO: Add VCF header according to VCF 4.2 standard
#TODO: Add INFO, FORMAT and SAMPLE column and its tags (values) according to 4.2 standard
#TODO: Add support for taking into account LIST of VCFs with known variants
#TODO: Quality of the base string and ref. genome will be added in later versions
