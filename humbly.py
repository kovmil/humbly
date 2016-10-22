#!/usr/bin/python
#------------------------------------------------------------------------------
# Filename: humbly.py
# Created: 06.10.2016
# Author: Milan Kovacevic, kovmil@gmail.com
#------------------------------------------------------------------------------
"""
Usage: humbly.py [-h] [--ctl=thr] [--clh=thr] [--known=vcf]... [--kvs=known_variants_significance] PILEUP_FILE

Options:
    -h --help
    --known Known variants file in uncompressed VCF format.
    --kvs Known variants significance constant [default: 1.2].
    --ctl Calling treshold low [default: 0.25]
    --cth Calling treshold high [default: 0.75]
"""
from docopt import docopt
import sys
import re
import mmap
from collections import Counter
from collections import namedtuple

SIGNIFICANT_PERCENTAGE_OF_BASES = (0.25)
MAJOR_PERCENTAGE_OF_BASES = (0.75)

FILTER_QUALITY = 50
FILTER_COVERAGE = 5
KNOWN_CONST = 1.2

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
    if len(alt) > 1:
        alt = alt[0]
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
    vcf += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    vcf += "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"
    vcf += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"

    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + str(sys.argv[1]).replace(".pileup", ".bam")

#------------------------------- Main program ---------------------------------


arguments = docopt(__doc__)
print(arguments)
#exit()

with open(arguments['PILEUP_FILE']) as pileup_file:
    lines = pileup_file.readlines()

#Variable for printing VCF file
header()

#C-like Typedef struct PiledupStruct
PileupStruct = namedtuple("PileupStruct", "chrom position ref coverage bases quality")

reg_exp = re.compile(r'([\.\^\$]*((\+|-)[0-9]+[ATCGN]+)*[ATCGN]+[\.\^\$]*)+')

for iter in lines:
    piled_up = str_to_pileup_struct(iter)

    mq_pos = str(piled_up.bases).find('^')
    if mq_pos >= 0:
        #ASCII value of next charachter - 33
        mapping_quality = ord(piled_up.bases[mq_pos+1]) - 33

    alt = None
    genotype = "0/0"
    Y_instead_indel = None

    #Regex for finding bases which are considered
    match_found = reg_exp.match(str(piled_up.bases))

    if match_found != None:
        base_len = len(str(match_found.group()))

        if  match_found.group().find('+') >= 0 or match_found.group().find('-') >= 0:
            if(match_found.group()[-1]!='.'):
                match = str(match_found.group())+"."
            else:
                match = str(match_found.group())

            if match.find('+') >= 0:
                first_plus_pos = match.find('+')
                num_of_indels = int(match[first_plus_pos+1])
                indel = match[first_plus_pos:first_plus_pos+num_of_indels+3]
            elif match.find('-') >= 0:
                first_minus_pos = match.find('-')
                num_of_indels = int(match[first_minus_pos+1])
                indel = match[first_minus_pos:first_minus_pos+num_of_indels+3]

            base_len = len(match.replace(indel, "Y"))
            Y_instead_indel = match.replace(indel, "Y")
            base_count =  Counter(match.replace(indel, "Y"))
        else:
            base_count = Counter(str(match_found.group()).replace("$","").replace("^",""))
            Y_instead_indel = str(piled_up.bases)


        mc_list = base_count.most_common(2)

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
        q = quality(Y_instead_indel.replace("$","").replace("^","") ,str(piled_up.quality), alt)

        known_snp = open('example_human_known_snp.vcf')
        known_snp_read = mmap.mmap(known_snp.fileno(), 0, access=mmap.ACCESS_READ)
        if known_snp_read.find("\t"+str(piled_up.position)+"\t") != -1:
            q = q * KNOWN_CONST

        indel_flag = None
        if alt == 'Y':
            indel_flag = "INDEL"
            alt = indel[:-1]
            known_indels = open('example_human_known_indels.vcf')
            known_indel_read = mmap.mmap(known_indels.fileno(), 0, access=mmap.ACCESS_READ)
            if known_indel_read.find("\t"+str(piled_up.position)+"\t") != -1:
                q = q * KNOWN_CONST

        #Printing VCF rows after header
        #Chrom
        vcf += "\n"+str(piled_up.chrom)
        #Position
        vcf += "  "+str(piled_up.position)
        #ID
        vcf += "  "+str(".")
        #Reference
        vcf += "  "+str(piled_up.ref)
        #Alt
        vcf += "  "+str(alt)
        #Quality
        vcf += "  "+str(q)
        #Filter
        if q > FILTER_QUALITY and int(piled_up.coverage) >= FILTER_COVERAGE:
            vcf += "  "+str("PASS")
        else:
            vcf += "  "+str(".")
        #Info
        vcf += "  "
        if indel_flag != None:
            vcf += "INDEL;"
        vcf += "DP="+str(piled_up.coverage)+";"
        vcf += "MQ="+str(mapping_quality)+";"
        #Format
        vcf += "  GT:PL"
        vcf += "  "+str(genotype)

print vcf
#TODO: Add VCF header according to VCF 4.2 standard
#TODO: Add INFO, FORMAT and SAMPLE column and its tags (values) according to 4.2 standard
