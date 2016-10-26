humbly - Simple(humble) Variant Caller

Python application used for calling variants on the aligned reads from the sequenced DNA.
It is constructred primarily for detecting, learning about bioinformatics and solving the challanges present in this process.

-Humbly can call both single nucleotide variants and indels 
-It can accept database of known variants 
-Configuration for various thresholds and filters supported 
-It is wrapped on Seven Bridges platform and tested against the existing variant callers achieving comperable results.

Usage: humbly.py [-h] [--ctl=thr] [--cth=thr] [--fq==filter_qual] [--fc==filer_coverage] [--known=vcf]... [--kvs=known_variants_significance] PILEUP_FILE
Options:
    -h --help
    --known Known variants file in uncompressed VCF format.
    --kvs Known variants significance constant [default: 1.2].
    --ctl Calling treshold low [default: 0.25]
    --cth Calling treshold high [default: 0.75]
    --fq Filter quality [default: 50]
    --fc Filter coverage [default: 5]
