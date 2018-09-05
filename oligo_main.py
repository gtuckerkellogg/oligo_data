import sys
sys.path.insert(0, "./oligo_py_import")

from oligo_get import *
from oligo_barcode import *


# TODO: use blockStarts in .bed file to display oligos for each gene
# TODO: implement amplification strategy-specific oligo design

# NOTE: working directory assumed to be /mnt/gtklab01/rachelle/oligo_data on atlas.cbis.nus.edu.sg


###################
##   VARIABLES   ##
###################

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

get_oligos(oligos = path_oligos,
           features = path_features,
           out_dir = path_out,
           name = project_name,
           strict = True, del_temp = False, temp_suff = "temp")


# NOTES ON EXECUTING 'assign_barcodes' FUNCTION
#   Assigns barcodes to oligos
# COMPULSORY INPUTS
#   oligos (str): path to .bed file of oligos (one oligo per line)
#   barcodes (str): path to file of barcodes
#   barcode_reader (func): function to read barcodes from file
#   out_dir (str): output directory
#   name (str): project name; to be used as name for output file
# OPTIONAL KEYWORD
#   append_barcode (bool): appends barcode to oligo if True, else prepends it
# OUTPUT
#   .bed file indexed by oligos (per file given to function) with barcode
#   id + sequence appended to each line of oligo

assign_barcodes(oligos = "EGR1_EGR1_oligos.bed",
                barcodes = "barcodes/bc25mer.240k.fasta",
                barcode_reader = read_DeLOB_barcodes,
                out_dir = "EGR1_barcodes_test",
                name = "EGR1_barcodes_test",
                append_barcode = True)










##############
##   TEST   ##
##############

test_all_features = {"A": "test_in/a.bed",
                 "B": "test_in/b.bed"}
test_selected_features = ["A", "B"]
test_features = [test_all_features[k] for k in test_all_features if \
                 k in test_selected_features]

test_oligos = "test_in/c.bed"
test_temp = "test-temp"
test_out = "test-out"
test_name = "test-strict-False"

# test
##get_oligos(oligos = test_oligos, features = test_features,
##           out_dir = test_out, name = test_name,
##           strict = False, del_temp = False, temp_suff = test_temp)




