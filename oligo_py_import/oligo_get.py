import shutil
import pybedtools

from oligo_errors import *


# get oligos with at least one nucleotide that overlaps with at least one feature in all
#   specified feature databases
# INPUTS
#   oligos (str): path to .bed file containing oligonucleotides
#   features (list): list of paths to .bed files containing features
#   out_dir (str): output directory
#   name (str): project name; to be used across all files/folders generated
# KEYWORDS
#   strict (bool): see docstring below for difference between True and False. Default value True
#   del_temp (bool): if True, temporary files created by function will be deleted at end of run. Defaut value is True
#   temp_suff (str): suffix for temporary folder. Default value is 'temp'.
# OUTPUTS
#   Writes file with name specified by variable 'name'. Contains sequences of
#   oligonucleotides with at least one nucleotide that overlaps with at least one
#   feature in all specified feature databases. Also contains genome coordinates of
#   oligonucleotides
def get_oligos(*, oligos, features, out_dir, name,
               strict = True, del_temp = True, temp_suff = "temp"):

    """
    Writes oligos with at least one nucleotide that overlaps with at least one feature
    in all specified feature databases to file.

    Note: If oligo_Z is 10 nt long and is checked against 2 feature files, A.bed and
    B.bed, and overlaps with feature A1 in feature file A.bed over the first
    3 nucleotides of oligo_Z, and with feature B1 in feature file B.bed over the last
    3 nucleotides of oligo_Z, it will NOT be reported if 'strict' is set to True.
    oligo_Z will be reported when 'strict' is set to True only if at least one nucleotide
    in oligo_Z has overlaps in ALL feature files. For more lenient selection, set 'strict'
    to False. (See below for difference)


    oligo_Z reported only when 'strict' is set to False

    oligo_Z                      ==========
    feature A1        ==============
    feature B1                          ==============
    shared overlaps   --------------------------------


    oligo_Z reported regardless of whether 'strict' is set to True or False.
        -Set 'strict' to True to report only oligos that fulfill this criterion

    oligo_Z                   ==========
    feature A1        ==============
    feature B1                    ==============
    shared overlaps --------------++------------
    """
    
    # read oligos .bed file into BedTool object
    oligos_bed = pybedtools.BedTool(oligos)
    # create output file path + name
    out_file = "{}{}{}.bed".format(out_dir, '/' if out_dir else '', name)
    
    # if out_path doesn't exist, create it
    if out_dir:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
    
    if strict:

        # PART 1: GET OLIGOS
        #   TODO: create function
        # create BedTool object from first feature database
        final_ranges = pybedtools.BedTool(features[0])
        
        # successive overlaps of databases to find regions that are shared
        #   across features from ALL databases
        for db in features[1:]:
            db_bed = pybedtools.BedTool(db)
            final_ranges = final_ranges.intersect(db_bed)

        # get appropriate oligos
        oligos_intersect = final_ranges.intersect(oligos_bed, wb = True)

        # create temporary directory to save intermediate results
        #   TODO: google temporary file creation in python (probably in sys)
        #       no need to worry about directory structure!
        temp_path = "{}_{}".format(name, temp_suff)
        if not os.path.exists(temp_path):
            os.makedirs(temp_path)

        # save intermediate results
        oligos_intersect.saveas("{}/{}_temp.bed".format(temp_path, name))


        # PART 2: REMOVE COORDINATES OF OVERLAPPED REGIONS (prettifying results)
        # open file of intersected oligos + create output file
        oligos_in = open("{}/{}_temp.bed".format(temp_path, name), 'r')
        oligos_out = open(out_file, "w+")

        # iterate through intersected oligos to remove coordinates of overlapped regions
        for entry in oligos_in.readlines():
            cols = entry.split('\t')[3:]
            oligos_out.write('\t'.join(cols))

        # close files
        oligos_in.close()
        oligos_out.close()

        # delete temporary directory if del_temp flag is raised
        if del_temp:
            shutil.rmtree("{}".format(temp_path))

    else:

        # successive overlaps of oligos with databases to find oligos that
        #   overlap with at least one feature in each database
        for db in features:
            db_bed = pybedtools.BedTool(db)
            oligos_bed = oligos_bed.intersect(db_bed, u = True)

        # save intermediate results
        oligos_bed.saveas(out_file)

    return
