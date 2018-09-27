import pybedtools

from oligo_gen import *


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
#   temp_dir (str): path to temporary directory (including name of temporary directory)
# OUTPUTS
#   Writes file with name specified by variable 'name'. Contains sequences of
#   oligonucleotides with at least one nucleotide that overlaps with at least one
#   feature in all specified feature databases. Also contains genome coordinates of
#   oligonucleotides
def get_oligos(*, oligos, features, out_dir, name,
               strict=True, del_temp=True, temp_dir=''):

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
    out_file = os.path.join(out_dir, name + ".bed")
    
    # if out_path doesn't exist, create it
    make_dir(dir_path=out_dir, check=False)
    
    if strict:

        # PART 1: GET OLIGOS
        
        # successive overlaps of databases to find regions that are shared
        #   across features from ALL databases
        common_ranges = get_common_ranges(*features)

        # get oligos that overlap with ranges in final_ranges
        oligos_intersected = common_ranges.intersect(oligos_bed, wb=True)

        # create temporary directory and save intermediate results
        temp_path = make_temp_dir(temp_dir=temp_dir,
                                  del_dir=del_temp,
                                  project_name=name)
        oligos_intersected.saveas(os.path.join(str(temp_path), name + ".bed"))

        # PART 2: REMOVE COORDINATES OF OVERLAPPED REGIONS (first 3 columns)
        remove_file_col(delimiter='\t',
                        fname_in=os.path.join(str(temp_path), name + ".bed"),
                        fname_out=out_file,
                        cols_to_remove=range(3))

        # CLEAN UP: delete temporary directory if del_temp flag is raised
        temp_path.clean_up()

    else:

        # successive overlaps of oligos with databases to find oligos that
        #   overlap with at least one feature in each database
        for db in features:
            db_bed = pybedtools.BedTool(db)
            oligos_bed = oligos_bed.intersect(db_bed, u=True)

        # save intermediate results
        oligos_bed.saveas(out_file)

    return


def get_common_ranges(*databases):

    """
    Accepts any number of paths to .bed files of database of features and outputs a
    BedTool object with ranges common to ALL databases given
    """
    
    if not databases:
        raise Exception("get_common_ranges requires at least 1 database.")

    common_ranges = pybedtools.BedTool(databases[0])

    # successive overlaps to find ranges common to all databases
    for db in databases[1:]:
        db_bed = pybedtools.BedTool(db)
        common_ranges = common_ranges.intersect(db_bed)

    return common_ranges
