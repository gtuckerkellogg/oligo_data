# from oligo_errors import *
from oligo_gen import *


###################
##     READ      ##
##   SEQUENCES   ##
###################


# all purpose function to read sequence files and store sequence as Seq objects in dict
#   indexed by sequence ids
def read_seq(fname, seq_coord, delimiter='\t', num_lines=1, id_gen=None) -> list:
    # open file containing sequences
    seqs = open(fname, 'r')

    # create output dict (seq_dict) and working variable (curr_seq)
    seq_out = []
    curr_seq = []

    # corral data into dictionary
    #   format: {<seq id>: [<seq>, <seq data>]}
    for i, line in enumerate(seqs.readlines()):

        # split current line by delimiter after removing line break character
        curr_seq.append(line[:-1].split(delimiter))

        # if lines belonging to current seq have all been read
        if len(curr_seq) == num_lines:

            # get seq id and sequence
            curr_id = i if not id_gen else id_gen(curr_seq)
            curr_sequence = curr_seq[seq_coord[0]][seq_coord[1]]

            # write seq to seq_dict
            seq_out.append(Seq(curr_id, curr_sequence, fname, curr_seq))
            curr_seq = []

    # close file
    seqs.close()

    return seq_out


def make_seq_reader(seq_coord, delimiter=None, num_lines=None, id_gen=None):
    def custom_seq_reader(fname):
        return read_seq(fname, seq_coord, delimiter=delimiter, num_lines=num_lines,
                        id_gen=id_gen)

    return custom_seq_reader


################
##   CUSTOM   ##
##   READERS  ##
################

# EXAMPLES!!!!!
# EXAMPLE FUNCTION TO READ SEQUENCES
# Note that both of the following two functions, 'read_seq_custom1' and 'read_seq_custom2',
# function the same way. Either way of generating a customised sequence reader is allowed.
def read_seq_custom1(fname):
    """
    Reads csv sequence file format where each line contains data for one sequence.
    Example of first 3 lines, where format is <id>,<seq>,<Tm> (Note, Tm values are bogus):
    id_1,ATGC,3.3
    id_2,AGGG,4.4
    id_3,GGCC,6.1
    """
    return read_seq(fname, [0, 1],
                    delimiter=',', id_gen=lambda x: x[0][0])


read_seq_custom2 = make_seq_reader([0, 1], delimiter=',', id_gen=lambda x: x[0][0])


# SPECIFIC FUNCTIONS FOR SPECIFIC FORMATS OF BARCODE FILES

# read barcodes generated by DeLOB
def read_DeLOB_barcodes(fname):
    return read_seq(fname, [1, 0],
                    num_lines=2, id_gen=lambda x: x[0][0])


def read_bedfile(fname):
    return read_seq(fname, [0, 3],
                    id_gen=lambda x: join_ele(',', x[0][:3]))


readers = {"bed": make_seq_reader([0, 3], id_gen=lambda entry: join_ele(',', entry[0][:3])),
           "delob": make_seq_reader([1, 0], num_lines=2, id_gen=lambda entry: entry[0][0])}


##################
##   SEQUENCE   ##
##    CLASS     ##
##################

class Seq:

    def __init__(self, id, seq, src, raw_data):
        self.id = id
        self.seq = seq
        self.src = src
        self.data = raw_data


class SeqFile:

    def __init__(self, fname, seq_list=None, file_type=None, file_reader=None, use_only=None):

        self.name = fname
        self.use_only = use_only
        self.ftype = file_type if file_type else os.path.splitext(fname)[0][1:]

        self.file_reader = file_reader if file_reader \
            else None if self.ftype not in readers \
            else readers[self.ftype]

        self.seq_list = seq_list if seq_list else \
            None if not file_reader else file_reader(fname)

    def fname(self) -> str:
        return self.name

    def get_seqs(self) -> list:
        return self.seq_list

    def num_seqs(self) -> int:
        return len(self.seq_list)

    def seq_at(self, i) -> Seq:
        if self.use_only:
            return self.seq_list[self.use_only]
        return self.seq_list[i]

    # sets file_reader and regenerates self.seq_list
    def set_file_reader(self, file_reader) -> None:
        self.file_reader = file_reader
        self.seq_list = file_reader(self.fname())

    # regenerates self.seq_list, and returns self.seq_list
    def read_seqs(self) -> list:
        if not self.get_seqs():
            raise Exception("No sequence reader found. Please use function 'make_seq_reader' to generate a reader." +
                            "Then use SeqFile.set_seq_reader() to assign reader to ths SeqFile object." +
                            "Finally, execute SeqFile.read_seqs() again.")
        self.seq_list = self.file_reader(self.fname())
        return self.seq_list
