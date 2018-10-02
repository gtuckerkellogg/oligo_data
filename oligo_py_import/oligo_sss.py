import math

from oligo_gen import *


#   .sss: sequence section source (arbitrary)
#       - stores final sequences + section starts + section descriptions + section ids + section sources
#       - Format (tab delimited):
#           - Col 1: main sequence
#           - Col 2: section 1 start position (wrt to <sequence, index from 0)
#           - Col 3: section 1 description
#           - Col 4: section 1 id (maybe not)
#           - Col 5: section 1 source file ID (arbitrarily assigned) (maybe not)
#           - Repeat col 2-5 for as many sub-sequences as there are in main sequence
#           - Each row contains a different main sequence
#           - Example:
#                   atgcatgcatgcatgc    0   hybridisation   oligo_01    source_001  4   long barcode    bcode_a11   source_002
#                   cgtacgtacgtacgta    0   hybridisation   oligo_02    source_001  4   long barcode    bcode_b389  source_002
#       - Last rows contain source information
#           - Each row corresponds with one source
#           - Col 1: source ID (arbitrarily assigned, but usually in order of appearance)
#           - Col 2: path to source file (not necessarily absolute)
#           - Example:
#                   source_001  C:/oligos.bed
#                   source_002  C:/barcodes.txt


################
##    FILE    ##
##   READER   ##
################

def read_sss(fname, subseq_cols=4):

    """
    Accepts a .sss file and reads it into an SSSFile object.
    """

    f = open(fname, 'r')

    seq_data = []
    src_data = {}

    for seq in f.readlines():

        cols = seq[:-1].split('\t')

        # if line describes source
        if cols[1][:6] == "source_":
            # src_data[cols[0]] = SSSSrc(cols[0], cols[1])
            src_data[cols[0]] = SSSSrc(cols[0], cols[1])

        # elif line describes sequence
        elif len(cols) > subseq_cols:
            seq_data.append({"main": cols[0],
                             "sub": [cols[i:i+subseq_cols] for i in range(1,
                                                                          len(cols)-subseq_cols,
                                                                          subseq_cols)]})

    f.close()

    return SSSFile(fname, raw_seq_data=seq_data, src_data=src_data)


#################
##   CLASSES   ##
#################

class SSSFile:

    def __init__(self, fname, raw_seq_data=None, sssseqs=None, src_data=None):

        self.name: str = fname
        self.sources_id: dict = {} if not src_data else src_data  # if provided, it's a dictionary indexed by source id (i.e. 'source_<some id>')
        self.sources_path: dict = invert_dict(self.sources_id)

        # store SSSSeq objects in list of provided, else start with empty list
        self.seqs: list = list(sssseqs) if sssseqs else []

        # if raw_seq_data instead of sssseqs, corral into SSSSeq objects
        if raw_seq_data and not sssseqs:
            for i, seq_entry in enumerate([] if not raw_seq_data else raw_seq_data):
                self.seqs.append(SSSSeq(seq_entry["main"],
                                        raw_data=seq_entry["subs"],
                                        order=i,
                                        sssfile=self))

    ''' GETTER FUNCTIONS (BASIC) '''

    def get_src_path(self, src_id) -> str:
        return None if (not self.sources_id and not self.sources_path) else\
            self.sources_id[src_id]

    def get_src_id(self, src_path) -> str:
        return None if (not self.sources_path and not self.sources_id) else\
            self.sources_path[src_path]

    ''' GETTER FUNCTIONS (PROCESSED) '''

    def search(self, query, *query_types, exact=False):

        """
        Accepts an arbitrary query and returns sequences with query phrase in any field
        (all fields are searched by default), but other 'query_types' options include
        "main", "sequence", "description", "source", and "id".
        Returns dictionary of relevant SSSSeq object indexed by row num in file.
        """
        
        query = str(query)
        query_types = list(query_types)
        output = {}
        query_types = ["main", "sequence", "description", "source", "id"] if not query_types \
                      else query_types

        # iterate through sequence entries
        for i, seq in enumerate(self.seqs):
            seq_data = []
            if "main" in query_types:
                query_types.remove("main")
                seq_data.append(seq.seq())
            seq_data.extend(list(map(lambda x: x.values(),
                                 seq.get_fields(*query_types).values())))
            if exact and query in seq_data:
                output[i] = seq
            elif not exact:
                for data_value in seq_data:
                    if query in data_value:
                        output[i] = seq
                        break
        return output

    def search_seq(self, query, sub_only=False, exact=False) -> dict:

        """
        Accepts a query for a sequence or main sequence, exact length or not.
        Returns dictionary of relevant SSSSeq object indexed by row num in file.

        To restrict search by excluding matches across different sub-sequences, set
        'sub_only' to True.
        To restrict search to exact matches (i.e. value in relevant field of positive
        results should have the exact same number of characters as the query in addition
        to containing the query), set 'exact' to True.
        """

        query = str(query)
        output = {}

        # iterate through sequence entries
        for i, seq in enumerate(self.seqs):
            if not sub_only:
                # search for matches in main sequence
                if (query == seq.seq() and exact) or \
                   (query in seq.seq() and not exact):
                        output[i] = seq
            else:
                # search for matches in sub-sequences
                for subseq in seq.subseqs:
                    if (query == subseq.seq() and exact) or \
                       (query in subseq.seq() and not exact):
                        output[i] = seq
                        break
        return output

    def seq_at_index(self, *indices, index_start=0) -> dict:

        """
        Get sequences at given index (i.e. row number). Returns a dictionary where keys
        are indices and values are SSSSeq objects.
        """

        output = {}
        for i in indices:
            i -= index_start
            try:
                output[i] = self.seqs[i]
            except IndexError:
                continue
        return output

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def add_seq(self, sssseq_obj) -> None:
        self.seqs.append(sssseq_obj)
        sssseq_obj.assign_sssfile(self)

    def recompile_sources(self) -> None:
        new_sources = [None]
        leading_zeroes = math.ceil(math.log(len(self.seqs), 10))

        # for each sequence in file
        for seq in self.seqs:
            # for each subsequence in file
            for subseq in seq.subseqs:
                # if path of source of current subsequence not already in new_sources
                if not subseq.src_path() in new_sources:
                    # add to path
                    new_sources.append(subseq.src_path())
                # generate new id name
                new_src_id = "source_{}".format(str(new_sources.index(subseq.src_path())).zfill(leading_zeroes))
                subseq.subseq_src_id = new_src_id
        self.sources_id = dict(("source_{}".format(str(i).zfill(leading_zeroes)), v)
                               for i, v in enumerate(new_sources))
        self.sources_path = invert_dict(self.sources_id)

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def write(self, fname, out_dir='') -> None:

        """
        Writes SSSFile object to .sss file.
        """

        make_dir(dir_path=out_dir, check=False)
        f = open(os.path.join(out_dir, fname), "w+")

        self.recompile_sources()

        for sequence in self.seqs:
            f.write(join_ele('\t', sequence.seq(), join_ele('\t',
                                                            *sequence.flatten_subseqs(source_id=True))) + '\n')

        for source in sorted([] if not self.sources_id else self.sources_id.items()):
            f.write(join_ele('\t', *source) + '\n')

        f.close()
        

class SSSSeq:

    def __init__(self, main, raw_data=None, ssssubseqs=None, sssfile=None, order=None):
        """
        :param main: main sequence composed of subsequences
        :param raw_data: raw data of subsequences, read almost directly from file (generally some iterable)
        :param ssssubseqs: subsequence data corraled into SSSSubSeq objects (also an iterable). Mutually exclusive with raw_data
        :param sssfile: file to which this SSSSeq object belongs
        :param order: order of SSSSeq in its SSSFile (i.e. row number)

        :type main: str
        :type sssfile: SSSFile
        :type order: int
        """

        self.seq_sequence: str = main
        self.seq_order: int = order
        self.sssfile: SSSFile = sssfile

        # store SSSSubSeq objects in list of provided, else start with empty list
        self.subseqs: list = list(ssssubseqs) if ssssubseqs else []   # list

        # if raw_data instead of ssssubseqs, corral into SSSSubSeq objects
        if not ssssubseqs:
            for i, sub_entry in enumerate(raw_data):
                # get start and end (given by start position of next sub-sequence)
                start = int(sub_entry[0])
                end = len(self.seq()) if (i >= len(raw_data - 1)) else int(raw_data[i+1][0])
                # create SubSeq object for current sub-sequence and add to output list
                self.subseqs.append(SSSSubSeq(self.seq_sequence[start:end], *sub_entry[:-1],
                                              src_id=sub_entry[-1], order=i))

        self.assign_sssseq_to_subseqs()

    def __repr__(self):
        return self.seq()

    ''' GETTER FUNCTIONS (BASIC) '''
    
    def seq(self) -> str:
        return self.seq_sequence

    def get_sssfile(self) -> SSSFile:
        return self.sssfile

    ''' GETTER FUNCTIONS (PROCESSED) '''

    def get_fields(self, *fields):

        """
        Receives any combination of valid field strings ("sequence", "start", "end",
        "description", "id", "source") and returns values of the specified fields
        (indexed by field name) of each subsequence (indexed by sub-sequence order).

        NOTE: this function takes advantage of the order-preserving property of lists and
        tuples. To be revised if this property changes.
        """

        field_types = {"sequence": lambda subseq: subseq.seq(),
                       "start": lambda subseq: subseq.start(),
                       "end": lambda subseq: subseq.end(),
                       "description": lambda subseq: subseq.description(),
                       "id": lambda subseq: subseq.id(),
                       "source": lambda subseq: subseq.src()}

        output = {}

        # iterate through all subsequences
        for subseq in self.subseqs:

            subseq_data = {}

            # get data for each field requested
            for field in fields:
                subseq_data[field] = field_types[field](subseq)

            # write data to output dictionary
            output[subseq.order] = subseq_data

        return output

    def get_subseqs_by_index(self, index_start=0, *i):
        """
        Accepts any number of integer values that are valid indices, returns a list of
        SSSSubSeq objects that are at the specified indices.
        :param i: valid indices, any quantity, integer values
        :param index_start: number to start counting indices from
        :return: list of SSSSubSeq objects at the specified indices (i)
        """
        output = []
        for subseq_order in i:
            output.extend(list(filter(lambda x: x.order() == subseq_order-index_start,
                                      self.subseqs)))
        return output

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def assign_sssfile(self, sssfile) -> None:
        if not self.sssfile:
            self.sssfile = sssfile
        else:
            raise Exception("This SSSSeq object already belongs to {}.".format(self.sssfile.name))

    def assign_sssseq_to_subseqs(self):
        for subseq in self.subseqs:
            subseq.assign_sssseq(self)

    def add_subseq(self, ssssubseq):
        """
        :type ssssubseq: SSSSubSeq
        """
        self.subseqs.append(ssssubseq)

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def flatten_subseqs(self, source_id=False):
        
        """
        Flattens sub-sequences data into single level list, returns the list. Each element is a new subsequence.
        Generally for use in writing to file.
        """

        output = []
        for subseq in self.subseqs:
            output.append(subseq.flatten(source_id=source_id))
        return output


class SSSSubSeq:

    # order starts from 0
    def __init__(self, seq, order, start, description, seq_id, src_id="source_None", src_path=None):

        self.subseq_sequence: str = seq
        self.subseq_order: int = order
        self.subseq_start: int = start
        self.subseq_description: str = description
        self.subseq_id: str = seq_id
        self.subseq_src_id: str = src_id
        self.subseq_src_path = src_path
        self.sssseq: SSSSeq = None

    def __repr__(self):
        return self.seq()

    ''' GETTER FUNCTIONS (BASIC) '''

    def id(self):
        return self.subseq_id

    def seq(self):
        """
        Returns string of self's sequence
        """
        return self.subseq_sequence

    def order(self):
        return self.subseq_order

    def start(self):
        return self.subseq_start

    def end(self):
        return self.subseq_start + len(self.seq()) - 1

    def pos(self):
        return self.start(), self.end()

    def description(self):
        return self.subseq_description

    def src_id(self):
        return self.subseq_src_id if self.subseq_src_id else \
            self.get_sssseq().get_sssfile().get_src_id(self.src_path()) if self.subseq_src_path else None

    def src_path(self):
        return self.subseq_src_path if self.subseq_src_path else \
            self.get_sssseq().get_sssfile().get_src_path(self.src_id()) if self.subseq_src_id else None

    def get_sssseq(self):
        return self.sssseq

    def get_sssfile(self):
        return self.get_sssseq().get_sssfile()

    ''' GETTER FUNCTIONS (PROCESSED) '''

    # 'source_id' indicates whether to substitute out path (str) with assigned ID of source
    def get_sssfile_fields(self, source_id=False):
        return self.start(), self.description(), self.id(), (self.src_id() if source_id else self.src_path())

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def assign_sssseq(self, new_sssseq):
        if not self.sssseq:
            self.sssseq = new_sssseq
        else:
            raise Exception("This subsequence already belongs to {}".format(self.sssseq))

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def flatten(self, source_id=False):
        return join_ele('\t', *self.get_sssfile_fields(source_id=source_id))


# unused
class SSSSrc:

    def __init__(self, src_id, src_path):
        self.src_id = src_id
        self.src_path = src_path

    def id(self):
        return self.src_id

    def path(self):
        return self.src_path
