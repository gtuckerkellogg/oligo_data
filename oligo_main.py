import sys
sys.path.insert(0, "./oligo_py_import")
sys.path.insert(0, "./amplification_strategies")
sys.path.insert(0, "./experiment_types")

from oligo_get import *
from oligo_design import *


# TODO: use blockStarts in .bed file to display oligos for each gene
# TODO: implement amplification strategy-specific oligo design
# TODO: create function to extract oligo sequences from .sss files

# NOTE: working directory assumed to be /mnt/gtklab01/rachelle/oligo_data on atlas.cbis.nus.edu.sg


###################
##   VARIABLES   ##
###################

# # # # Variables for generating file of TARGET HYBRIDISATION SEQUENCES
# # Path variables should be modified according to working directory
# dictionary of paths to feature databases
all_features = {"EGR1": "../OligoMiner/20180228_EGR1_oligominer_u/egr1_promoter.bed",
                "karpas_DMSO_H3K4me3": "../karpas_write/MACS/H3K4me3/20180813_DMSO/20180813_DMSO_peaks.bed",
                "karpas_DMSO_H3K27me3": "../karpas_write/MACS/H3K27me3/20171205_DMSO_1/20171205_DMSO_1_peaks.bed"}
# select databases to use (this is a list of keys of databases in all_features dict)
selected_features = ["EGR1", "EGR1"]

# other modifiable variables
path_oligos = "/mnt/gtklab01/shared/oligominer/hg38b/hg38_chr5b.bed"
path_temp = "temp"  # don't leave this empty
path_out = ''       # leave as '' if output file is to be written to working directory, else provide path to desired output directory
project_name = "EGR1_EGR1_oligos"

# DO NOT MODIFY the following line
path_features = [all_features[k] for k in all_features if k in selected_features]


# # # # Variables for AMPLIFICATION STRATEGIES
# can be used to build new strategy templates
amp_templates = {"PEA_PCR": "amplification_strategies/PEA_PCR.txt",
                 "PLA_PCR": "amplification_strategies/PLA_PCR.txt",
                 "RC_PLA": "amplification_strategies/RC_PLA.txt"}

# used to generate oligonucleotides
amp_strategies = {"PEA_PCR": "amplification_strategies/PEA_PCR_example.txt",
                  "PLA_PCR": "amplification_strategies/PLA_PCR_example.txt",
                  "RC_PLA": "amplification_strategies/RC_PLA_example.txt",
                  "example": "amplification_strategies/example.txt"}


####################
##   EXECUTIONS   ##
####################

# NOTES ON EXECUTING 'get_oligos' FUNCTION
# COMPULSORY FIELDS
#   oligos (str): path to .bed file containing oligonucleotides
#   features (list): list of paths to .bed files containing features
#   out_dir (str): output directory
#   name (str): project name; to be used across all files/folders generated
# KEYWORDS
#   strict: default value True. See docstring difference between True and False
#   del_temp: default value True. If True, temporary files created by function will be deleted at end of run
#   temp_suff: default value 'temp'. Suffix for temporary folder.
# OUTPUTS
#   Writes file with name specified by variable 'name'. Contains sequences of
#   oligonucleotides with at least one nucleotide that overlaps with at least one
#   feature in all specified feature databases. Oligonucleotide sequences are indexed
#   by genomic coordinates

# get oligos that overlap with features of interest
get_oligos(oligos=path_oligos,
           features=path_features,
           out_dir=path_out,
           name=project_name,
           strict=True, del_temp=False, temp_dir=path_temp)


# if FROM_READYMADE is set to True, RC_PLA_A will be obtained from ready-made
#   amplification strategy from 'amp_strategies'. Else, RC_PLA_A will be
#   built from empty template. Either way, the final RC_PLA_A design is identical
FROM_READYMADE = False

# starting from ready-made oligo design
if FROM_READYMADE:

    # read oligo design for desired strategy into SODFile object
    RC_PLA = read_sod(amp_strategies["RC_PLA"])
    # extract out SODOligo object containing design for desired oligo (in this case, oligo A)
    RC_PLA_A = RC_PLA.get_oligo('A')

# creating new oligo design based on empty template
else:

    # read empty oligo design template for desired strategy into SODFile object
    RC_PLA = read_sod(amp_templates["RC_PLA"])
    # extract out SODOligo object containing design for desired oligo (in this case, oligo A)
    RC_PLA_A = RC_PLA.get_oligo('A')

    # provide sequences for blank fields
    RC_PLA_A.replace_seq("A_3_1", "LINKERSEQUENCE")
    RC_PLA_A.replace_seq("A_3_2", "CIRCLE_PART1_5'_COMP")
    RC_PLA_A.replace_seq("A_3_4", "CIRCLE_PART2_3'_COMP")

    
# obtains blank sub-sequences and prints blank sub-sequences to be filled in
#   (in this context, this statement is for viewing only and doesn't perform
#   any useful function)
RC_PLA_A_blank = RC_PLA_A.blank_subseqs(show=True)

# create OligoBuilder object by feeding it the SODOligo object for the desired oligo
RC_PLA_A_builder = OligoBuilder(RC_PLA_A)

# read files with sequences to use into SeqFile objects
#   Note that file readers provided in this example are defined in oligo_seq.py
#   Instructions for customising file readers can be found in oligo_seq.py
EGR1_hyb_seq = SeqFile("EGR1_EGR1_oligos.bed", file_reader=read_bedfile)
barcodes = SeqFile("barcodes/bc25mer.240k.fasta", file_reader=read_DeLOB_barcodes)

# specify which SeqFile object is for which sub-sequence
RC_PLA_A_builder.add_seqfile("A_0_0", EGR1_hyb_seq)
RC_PLA_A_builder.add_seqfile("A_3_3", barcodes)

# create oligos and write to .sss file, TADAA!
RC_PLA_A_SSS = RC_PLA_A_builder.create_oligos("EGR1_barcodes_RC_PLA_A.sss", out_dir="test_out")

# RC_PLA_A_builder can be reused by providing different SeqFile objects
# by re-executing OligoBuilder.add_seqfile() and OligoBuilder.create_oligos()
