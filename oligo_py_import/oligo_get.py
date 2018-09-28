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
               strict=True, del_temp=True, temp_dir='', sorted=True):

    """
    Writes oligos with at least one nucleotide that overlaps with at least one feature
    in all specified feature databases to file.

    Assumes that files are sorted, which invokes a more efficient algorithm. If entries
    in files are not sorted, set 'sorted' to False.

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
        
        # successive overlaps of databases to find regions that are shared
        #   across features from ALL databases
        common_ranges = get_common_ranges(*features)

        # get oligos that overlap with ranges in final_ranges
        oligos_bed = oligos_bed.intersect(common_ranges, u=True, sorted=sorted)

    else:

        # successive overlaps of oligos with databases to find oligos that
        #   overlap with at least one feature in each database
        for db in features:
            db_bed = pybedtools.BedTool(db)
            oligos_bed = oligos_bed.intersect(db_bed, u=True)

    # save intermediate results
    oligos_bed.saveas(os.path.join(out_dir, out_file))

    return


def pos_to_block_start(feature_data, *pos_coord) -> str:
    """
    Accepts any quantity of (<start>, <end>) tuple pairs that describe a block within
    a bedfile feature.
    :param feature_data: dictionary of feature data from bedfile
    :param pos_coord: (<start>, <end>) tuple pairs, any quantity
    :return: formatted string for bedfile entry of feature with specified blockstarts
    """

    # corral feature_data into dictionary to enable use of dict.get
    feature_dict = dict((i, v) for i, v in enumerate(feature_data))

    # get start and end positions of all blocks
    starts, ends = zip(*pos_coord)

    # provide default values for bedfile fields
    fields_1 = ["Unknown", int(starts[0]), int(ends[-1]), "Unknown", 0, '+']
    fields_2 = (0, 1,
                max(int(feature_dict[2]), int(ends[-1])) - min(int(feature_dict[1]), int(starts[0])) + 1,
                0)

    # create bedfile entry (str)
    return join_ele('\t', *map(lambda i: feature_dict.get(i, fields_1[i]), range(len(fields_1))),
                    join_ele(',', *starts),
                    join_ele(',', *ends),
                    *fields_2)


def make_oligo_blocks(fname, features, oligos, sorted=True):

    features_bed = pybedtools.BedTool(features)
    oligos_bed = pybedtools.BedTool(oligos)

    out_bed = []

    for feature in features_bed:

        # create temporary bedtool object using feature to facilitate intersection, then intersect
        feature_bed = pybedtools.BedTool(join_ele('\t', *feature), from_string=True)
        feature_oligos = oligos_bed.intersect(feature_bed, u=True, sorted=sorted)

        # retrieve coordinates of oligos
        feature_oligos_coords = [tuple(int(i) for i in oligo[1:3]) for oligo in feature_oligos]

        # generate bedfile entry
        out_bed.append(pos_to_block_start(list(feature), *feature_oligos_coords))

    out_bedfile = pybedtools.BedTool(join_ele('\n', out_bed), from_string=True)
    out_bedfile.saveas(fname)


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
