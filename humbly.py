#!/usr/bin/python
#------------------------------------------------------------------------------
# Filename: humbly.py
# Created: 06.10.2016
# Author: Milan Kovacevic, kovmil@gmail.com
#------------------------------------------------------------------------------
"""
Usage: humbly.py [-h] [--ctl=thr] [--cth=thr] [--fq==filter_qual] [--fc==filer_coverage] [--known=vcf]... [--kvs=known_variants_significance] PILEUP_FILE

Options:
    -h --help
    --known Known variants file in uncompressed VCF format.
    --kvs Known variants significance constant [default: 1.2].
    --ctl Calling treshold low [default: 0.25]
    --cth Calling treshold high [default: 0.75]
    --fq Filter quality [default: 50]
    --fc Filter coverage [default: 5]
"""
from docopt import docopt
import sys
import re
import mmap
from collections import Counter
from collections import namedtuple

FILTER_QUALITY = 50
FILTER_COVERAGE = 5
KNOWN_CONST = 1.2
CALLING_TRESHOLD_HIGH = 0.8
CALLING_TRESHOLD_LOW = 0.2

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
    # print base
    # print qual
    # print alt
    # print "-----------------------------"
    if len(alt) > 1:
        alt = alt[0]
    if abs(len(base) - len(qual)) == 1:
        base = base[:-1]
    if len(base) == len(qual):
        positions = [m.start() for m in re.finditer(str(alt), base)]
        num_of_chosen = 1
        if len(positions) > 0:
            num_of_chosen = len(positions)
            s = 0
            for n in positions:
                s += ord(qual[n]) #ord(c) returns ASCII value of 'c'
            return round(s*1.0/num_of_chosen, 2)
    else:
        return 0

def header():
    global vcf
    vcf += "##fileformat=VCFv4.2\n"
    vcf += "##FILTER=<ID=PASS, Description=\"All filters passed\">\n"
    vcf += "##ALT=<ID=*,Description=\"Represents allele(s) other than observed.\">\n"
    vcf += "##INFO=<ID=INDEL,Number=0,Type=Flag,Description=\"Indicates that the variant is an INDEL.\">\n"
    vcf += "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Raw read depth\">\n"
    vcf += "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"List of Phred-scaled genotype likelihoods\">\n"
    vcf += "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"

    vcf += "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + str(arguments['PILEUP_FILE']).replace(".pileup", ".bam")

#------------------------------- Main program ---------------------------------

arguments = docopt(__doc__)
#print(arguments)

if arguments['--cth'] == None:
    call_thr_high = CALLING_TRESHOLD_HIGH
else:
    call_thr_high = int(arguments['--cth'])
if arguments['--ctl'] == None:
    call_thr_low = CALLING_TRESHOLD_LOW
else:
    CALLING_TRESHOLD_LOW = int(arguments['--ctl'])

if arguments['--fq'] == None:
    filter_qual = FILTER_QUALITY
else:
    filter_qual = int(arguments['--fq'])

if arguments['--fc'] == None:
    filter_cov = FILTER_COVERAGE
else:
    filter_cov = int(arguments['--fc'])

with open(arguments['PILEUP_FILE']) as pileup_file:
    lines = pileup_file.readlines()

#Function for creating VCF header
header()

# C-like Typedef struct PiledupStruct
PileupStruct = namedtuple("PileupStruct", "chrom position ref coverage bases quality")

# Compile regular expression for finding nucleotides in pileup
reg_exp = re.compile(r'.*[ATCGN]+.*')


if arguments['--known'] != []:
    known_snp = open(str(arguments['--known'][0]))
    known_snp_read = mmap.mmap(known_snp.fileno(), 0, access=mmap.ACCESS_READ)

for iter in lines:
    piled_up = str_to_pileup_struct(iter) # Pileup line to structure


    mq_pos = str(piled_up.bases).find('^')
    if mq_pos >= 0:
        #ASCII value of next charachter - 33
        mapping_quality = ord(piled_up.bases[mq_pos+1]) - 33

    alt = None
    genotype = "0/0"
    Y_instead_indel = None

    #Regex for finding bases which are considered
    match_found = reg_exp.match((str(piled_up.bases).replace("$","").replace("^","")))

    if match_found != None:
        base_len = len(str(match_found.group()))

        # Prepare for detecting indels
        indel_flag = None
        if  match_found.group().find('+') >= 0 or match_found.group().find('-') >= 0:
            indel_flag = "INDEL"
            if(match_found.group()[-1]!='.'):
                match = str(match_found.group())+"."
            else:
                match = str(match_found.group())

            if match.find('+') >= 0 and piled_up.bases[piled_up.bases.find('+')+1].isdigit():
                first_plus_pos = match.find('+')
                num_of_indels = int(match[first_plus_pos+1])
                indel = match[first_plus_pos:first_plus_pos+num_of_indels+3]
            elif match.find('-') >= 0 and piled_up.bases[piled_up.bases.find('-')+1].isdigit():
                first_minus_pos = match.find('-')
                num_of_indels = int(match[first_minus_pos+1])
                indel = match[first_minus_pos:first_minus_pos+num_of_indels+3]
            else:
                Y_instead_indel = match.replace("+","").replace("-","")

            Y_instead_indel = match.replace(indel, "Y")

            while Y_instead_indel.find('+')>=0 or Y_instead_indel.find('-')>=0:
                indel_flag = None
                if Y_instead_indel.find('+')>=0 and piled_up.bases[piled_up.bases.find('+')+1].isdigit():
                    first_plus_pos = Y_instead_indel.find('+')
                    num_of_indels = int(Y_instead_indel[first_plus_pos+1])
                    indel = Y_instead_indel[first_plus_pos:first_plus_pos+num_of_indels+3]
                    Y_instead_indel = Y_instead_indel.replace(indel, "Z")
                elif Y_instead_indel.find('-') >= 0 and piled_up.bases[piled_up.bases.find('-')+1].isdigit():
                    first_minus_pos = Y_instead_indel.find('-')
                    num_of_indels = int(Y_instead_indel[first_minus_pos+1])
                    indel = Y_instead_indel[first_minus_pos:first_minus_pos+num_of_indels+3]
                    Y_instead_indel = Y_instead_indel.replace(indel, "Z")
                else:
                    Y_instead_indel = match.replace("+","").replace("-","")

            base_len = len(Y_instead_indel)

            base_count =  Counter(Y_instead_indel.replace("$","").replace("^",""))
        else:
            base_count = Counter(str(match_found.group()).replace("$","").replace("^",""))
            Y_instead_indel = str(piled_up.bases)

        # Finding most common and second most common variants in the pileup
        mc_list = base_count.most_common(2)

        mc = mc_list[0][0]
        mc_cnt = mc_list[0][1]

        if len(mc_list) > 1:
            smc = mc_list[1][0]
            smc_cnt = mc_list[1][1]
        else:
            smc = None
            smc_cnt = 0

        known_found = -1
        posit = str(piled_up.position)
        # print posit
        # print smc
        if arguments['--known'] != []:
            known_found = known_snp_read.find("\t"+posit+"\t")


        if mc == '.':
            if smc_cnt > base_len * call_thr_low and smc_cnt>1:
                alt = smc
                genotype = "0/1"
            else:
                # if known_found != -1 and smc_cnt>base_len * 0.2:
                #     alt = smc
                #     genotype = "0/1"
                # else:
                alt = piled_up.ref
                genotype = "0/0"
        else:
            alt = mc
            if mc_cnt > base_len * call_thr_high and mc_cnt>1:
                genotype = "1/1"
            elif smc_cnt > base_len * call_thr_low and mc_cnt>1:
                if smc == '.':
                    genotype = "0/1"
                else:
                    alt = [alt, smc]
                    genotype = "1/2"
            else:
                alt = piled_up.ref
                genotype = "0/0"


    if genotype != "0/0":
        variant_quality = quality(Y_instead_indel.replace("$","").replace("^","").replace("]", "").replace("!", "").replace("I", "") ,str(piled_up.quality), alt)

        indel_flag = None

        # if arguments['--known'] != []:
        #     if known_snp_read.find("\t"+str(piled_up.position)+"\t") != -1:
        #         variant_quality = variant_quality * KNOWN_CONST


        reference = piled_up.ref
        if alt == 'Y' or alt == 'Z':
            indel_flag = "INDEL"
            if indel[0] == '+':
                alt = reference + indel[2:-1]
            else:
                alt = reference
                reference = reference + indel[2:-1]

        #Printing VCF rows after header
        #Chrom
        vcf += "\n"+str(piled_up.chrom)
        #Position
        vcf += "\t"+str(piled_up.position)
        #ID
        vcf += "\t"+str(".")
        #Reference
        vcf += "\t"+str(reference)
        #Alt
        if genotype == "1/2":
            vcf += "\t"+str(alt[0])+","+str(alt[1])
        else:
            vcf += "\t"+str(alt)
        #Quality
        if variant_quality == None:
            variant_quality = 0
        vcf += "\t"+str(variant_quality)
        #Filter
        if variant_quality > filter_qual and int(piled_up.coverage) >= filter_cov:
            vcf += "\t"+str("PASS")
        else:
            vcf += "\t"+str(".")
        #Info
        vcf += "\t"
        if indel_flag != None:
            vcf += "INDEL;"
        vcf += "DP="+str(piled_up.coverage)+";"
        vcf += "MQ="+str(mapping_quality)+";"
        #Format
        vcf += "\tGT"
        vcf += "\t"+str(genotype)

#Printing the output VCF
output_file = open('output.vcf', 'w')
output_file.write(vcf)
output_file.close()

#DEBUG MODE
#print vcf

#TODO: Include quality into variant/genotype decisioning
#TODO: Add INFO, FORMAT and SAMPLE column and its tags (values) according to 4.2 standard
#TODO: output file name
