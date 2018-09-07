from oligo_errors import *
from oligo_gen import *

# TODO: create file format for barcoded oligos


#################
##   BARCODE   ##
##    CLASS    ##
#################

class Barcode():

    def __init__(self, id, seq, raw_data):
        self.id = id
        self.seq = seq
        self.data = raw_data


##################
##     READ     ##
##   BARCODES   ##
##################



# all purpose function to read barcode files and store barcodes as Barcode objects in dict
#   indexed by barcode ids
def read_barcodes(fname, delimiter, num_lines, barcode_coord, id_coord):

    # open file containing barcodes
    barcodes = open(fname, 'r')

    # create output dict (barcode_dict) and working variable (curr_barcode)
    barcode_dict = {}
    curr_barcode = []

    # corral data into dictionary
    #   format: {<barcode id>: [<barcode>, <barcode data>]}
    for line in barcodes.readlines():
        
        # split current line by delimiter after removing line break character
        curr_barcode.append(line[:-1].split(delimiter))

        # if lines belonging to current barcode have all been read
        if len(curr_barcode) == num_lines:
            
            # get barcode id and sequence
            curr_id = curr_barcode[id_coord[0]][id_coord[1]]
            curr_seq = curr_barcode[barcode_coord[0]][barcode_coord[1]]
            
            # write barcode to barcode_dict
            barcode_dict[curr_id] = Barcode(curr_id, curr_seq, curr_barcode)
            curr_barcode = []

    # close file
    barcodes.close()

    return barcode_dict




## EXAMPLE FUNCTION TO READ BARCODE
def read_barcode_custom1(fname):
    """
    Reads csv barcode file format where each line contains data for one barcode.
    Example of first 3 lines, where format is <id>,<seq>,<Tm> (Note, Tm values are bogus):
    id_1,ATGC,3.3
    id_2,AGGG,4.4
    id_3,GGCC,6.1
    """
    return read_barcode(fname, ',', 1, [0,1], [0,0])


## SPECIFIC FUNCTIONS FOR SPECIFIC FORMATS OF BARCODE FILES

# read barcodes generated by DeLOB
def read_DeLOB_barcodes(fname):
    return read_barcodes(fname, '\t', 2, [1,0], [0,0])




##################
##    ASSIGN    ##
##   BARCODES   ##
##################

# assigns barcodes to oligos
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
def assign_barcodes(*, oligos, barcodes, barcode_reader, out_dir, name,
                    append_barcode = True):

    # read barcodes, create output file path + name
    barcode_list = list(barcode_reader(barcodes).values())
    out_file = "{}{}{}.bed".format(out_dir, '/' if out_dir else '', name)
    
    # if out_path doesn't exist, create it
    if out_dir:
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)

    # open file of oligos + create output file
    oligos_data = open(oligos, 'r')
    oligos_out = open(out_file, "w+")

    i = 0

    # append barcode id + sequence to each oligo entry in output file
    for oligo in oligos_data.readlines():
        # TODO: make getting oligo sequence more robust
        oligo_seq = oligo.split('\t')[-2]
        barcode = barcode_list[i]
        final_seq = (oligo_seq + barcode.seq) if append_barcode else \
                    (barcode.seq + oligo_seq)
        oligos_out.write("{}\t{}\t{}\t{}\n".format(oligo[:-1], barcode.id,
                                                   barcode.seq, final_seq))
        i += 1

    # close files
    oligos_data.close()
    oligos_out.close()

    return
    
