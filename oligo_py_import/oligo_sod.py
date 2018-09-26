import re

from oligo_gen import *
from oligo_sss import *
from oligo_seq import *
 
# .sod: Sequences of Oligo Design
#   - stores template for constructing oligos
#   - Format (tab delimited):
#       - Header: ">> .SOD\n"
#          - Indicates that format is compatible with SODFile object. Required as first line.
#       - Col 1: sub-sequence id
#          - <oligo_id>_<pole>_<distance from 0>
#              - oligo_id: id of oligonucleotide (e.g. 'A')
#              - pole: indicates whether sub-seq is at 3' (3) or 5' (5) end or is central (0)
#              - distance from 0: indicates position of sub-seq relative to main (some integer)
#       - Col 2: sequence of sub-sequence, or <> or <<>> if currently undefined
#          - '<>' if field is compulsory for building oligo
#          - '<<>>' if field is optional for building oligo
#       - Col 3: description of sub-sequence
#       - Each row contains a different sub-sequence. Doesn't have to be in order.
#       - Lines of comments begin with '#'
#       - Example:
#               >> .SOD
#
#               A_5_2  B       Biotin
#               A_5_1  <<>>    linker sequence
#               A_0_0	<>	target hybridisation sequence
#               A_3_1	<<>>	linker sequence
#               A_3_2	<>	target barcode
#               A_3_3	<>	sequence complementary to 3' of oligo B
#
#               B_0_0	<<>>	target hybridisation sequence
#               B_3_1	<<>>	linker sequence
#               B_3_2	<>	target barcode or common primer
#               B_3_3	<>	sequence complementary to 3' of oligo A
#
#               # A: main oligo
#               # B: either attached to Ab or another target

# TODO: function to generate oligos and write to file (.sss?)


################
##    FILE    ##
##   READER   ##
################

def read_sod(fname):

    """
    Accepts a .sod file that contains information on subsequences of relevant
    oligonucleotides for a given experiment type/amplification strategy and
    reads it into a SODFile object.
    """
    f = open(fname, 'r')

    oligos = {}
    lines = f.readlines()

    if lines[0] != ">> .SOD\n":
        print(lines[0])
        raise Exception("Invalid file format detected. Please try another file.")

    for line in lines:

        # if not comments (i.e. if line doesn't start with '#')
        if line[0] != '#' and line[0] != '>':
            
            cols: list = line[:-1].split('\t')

            # if there are at least 2 columns (i.e. not blank line)
            if len(cols) >= 2:

                if len(cols) == 2:
                    cols.append('No description')
                
                # get the relevant id values + sequence/description info
                #   subseq_id is in <oligo id>_<pole (3' or 5' or 0)>_<distance from 0> format
                subseq_id = SODId(cols[0])
                subseq_obj = SODSubSeq(subseq_id, cols[1], join_ele('\t', *cols[2:]))
                subseq_id.set_subseq(subseq_obj)

                # write to oligos dictionary
                oligo_id = subseq_id.oligo_id
                oligos[oligo_id] = oligos.get(oligo_id, []) + [subseq_obj]
    
    # create SODFile ojbect and return
    return SODFile(fname, oligos)


#################
##   CLASSES   ##
#################

class SODFile:

    def __init__(self, fname, oligos_data):

        self.name = fname

        # corral oligos_data into dictionary index subsequences by oligo ids
        #   structure:  {'A': <SODOligo object 1>,
        #                'B': <SODOligo object 2>}
        self.oligos = {}
        for oligo_id, oligo_subseqs in oligos_data.items():
            self.oligos[oligo_id] = SODOligo(oligo_id, oligo_subseqs, sodfile=self)

    def fname(self) -> str:
        return self.name

    def oligo_ids(self) -> list:

        """
        Return list of oligo ids
        """
        return list(self.oligos.keys())

    def valid_oligo_id(self, oligo_id) -> bool:
        return oligo_id in self.oligo_ids()

    def num_oligos(self) -> int:
        return len(self.oligos)

    # retrieves relevant SODOligo object
    def get_oligo(self, oligo_id):
        """
        :param oligo_id: SODOligo ID (str)
        :return: SODOligo object
        """
        return self.oligos[oligo_id]
    
    def get_3prime_design(self, oligo_id) -> list:
        return self.oligos[oligo_id].get_3prime_design()

    def get_5prime_design(self, oligo_id) -> list:
        return self.oligos[oligo_id].get_5prime_design()

    def get_0(self, oligo_id) -> list:
        return self.oligos[oligo_id].get_0()

    def get_subseq(self, subseq_id):
        """
        :param subseq_id: SODSubSeq ID, as string ('A_0_0') or tuple (('A', 0, 0)) or SODId object
        :return: SODSubSeq object that has the given ID
        """
        oligo_obj = self.get_oligo(SODId(subseq_id).oligo_id)
        return oligo_obj.get_subseq(subseq_id)

    def replace_seq(self, subseq_id, subseq_seq) -> None:
        oligo_obj = self.get_oligo(SODId(subseq_id).oligo_id)
        oligo_obj.replace_seq(subseq_id, subseq_seq)

    # retrieves design of relevant oligo as list of SODSubSeq objects.
    def get_design(self, oligo_id) -> list:
        return self.oligos[oligo_id].get_design()

    def build_oligo(self, oligo_id, *seqs, sss=False) -> str or SSSSeq:
        return self.oligos[oligo_id].build_oligo(*seqs, sss=sss)

    def write(self, fname, out_dir='', comments='') -> None:

        """
        Writes SODFile object to .sod file.
        """

        make_dir(dir_path=out_dir, check=False)
        f = open(os.path.join(out_dir, fname), "w+")

        oligos = list(map(lambda x: x[1],
                          sorted(self.oligos.items(),
                                 key=lambda x: x[0])))

        # create ">> .SOD" header
        f.write(">> .SOD\n\n")

        # write sub-sequence data
        for i, oligo in enumerate(oligos):
            f.write(join_ele('\n', *list(map(lambda subseq: subseq.flatten(),
                                             oligo.get_design()))) +\
                    ('' if i >= len(oligos) else '\n\n'))

        # write comments
        if isinstance(comments, list) or isinstance(comments, tuple):
            comments = join_ele('\n', *comments)
        comments = replace_ele(comments, '\n', "\n# ")
        f.write("# " + comments)

        f.close()


class SODOligo:

    def __init__(self, oligo_id, subseqs, sodfile=None):

        self.oligo_id: str = oligo_id

        # list of SODSubSeq objects
        self.subseqs: list = subseqs
        # set parent oligos of subseq_id (SODId) and subseq_obj (SODSubSeq)
        for subseq_obj in self.subseqs:
            subseq_obj.set_oligo(self)
            subseq_obj.id().set_oligo(self)

        # SODFile object to which this SODOligo belongs
        self.sodfile: SODFile = sodfile

    def id(self) -> str:
        return self.oligo_id

    # Set SODFile object to which this SODOligo belongs
    def set_sodfile(self, sodfile) -> None:
        self.sodfile = sodfile

    def file(self) -> SODFile:
        return self.sodfile

    # prints sub-sequence fields in pretty-ish format
    #   supported fields: "id" ("i"), "sequence" ("seq", "s"), "description" ("desc", "d")
    def pretty_print_subseqs(self, *fields) -> None:
        
        header = []
        header += ["ID"] if contains_any(fields, "ID", "id", 'i') else []
        header += ["Sequence"] if contains_any(fields, "Sequence", "sequence", "seq", 's') else []
        header += ["Description"] if contains_any(fields, "Description", "description", "desc", 'd') else []
        print(join_ele('\t', *header))

        for subseq_obj in self.subseqs:
            subseq_obj.pretty_print(*header)

    def filter_id(self, pole: int = None, distance: int = None) -> list:

        """
        Returns list of SODSubSeq objects of subsequences at a given pole
        (5' (pole = 5) or 3' (pole = 3) or none (pole = 0)) and/or distance from 0.
        Both 'pole' and 'distance' variables only support integer values.
        Not sorted in any order.
        """
        output = []
        for subseq_obj in self.subseqs:
            if (pole is None or subseq_obj.pole == pole) and\
               (distance is None or subseq_obj.distance == distance):
                output.append(subseq_obj)
        return output

    def get_subseq(self, subseq_id):

        """
        Accepts subseq_id in 'A_0_0' or ('A', 0, 0) or SODId formats.
        :return: SODSubSeq object if found, else None
        """
        for subseq_obj in self.subseqs:
            if subseq_obj.id() == subseq_id:
                return subseq_obj
        print("{} not found in oligo {}.".format(str(subseq_id), self.id()))

    def get_3prime_design(self) -> list:
        
        """
        Returns list of SODSubSeq objects of subsequences at the 3' end.
        Not sorted in any order.
        """
        return self.filter_id(pole=3)

    def get_5prime_design(self) -> list:
        
        """
        Returns list of SODSubSeq objects of subsequences at the 5' end.
        Not sorted in any order.
        """
        return self.filter_id(pole=5)

    def get_0(self) -> list:
        
        """
        Returns list of SODSubSeq object of subsequence at 0.
        """
        return self.filter_id(pole=0)

    def get_design(self) -> list:

        """
        Returns a list of SODSubSeq objects of subsequences in 5' --> 3' order.
        """
        five_prime = sort_sodsubseqs(*self.get_5prime_design(), reverse=True)
        three_prime = sort_sodsubseqs(*self.get_3prime_design())
        main = self.get_0()
        return list(five_prime + main + three_prime)

    def blank_subseqs(self, show=False) -> list:

        """
        Returns a list of SODSubSeq objects that can be filled in upon execution of
        self.build_oligo()

        Set 'show' to True to print ID and description of blank sub-sequences.
        Use to identify IDs of blanks and their purpose, so that the appropriate
        sequences can be provided when calling self.build_oligo()
        
        TODO: maybe combine with build_oligos?
        """
        output = []
        for subseq_obj in self.subseqs:
            if subseq_obj.is_blank():
                if show:
                    subseq_obj.pretty_print("id", "description")
                output.append(subseq_obj)
        return output

    def compulsory_subseqs(self, show: bool = False) -> list:
        
        """
        Returns a list of SODSubSeq objects that MUST be filled in upon execution of
        self.build_oligo()

        Set 'show' to True to print ID and description of blank sub-sequences.
        Use to identify IDs of blanks and their purpose, so that the appropriate
        sequences can be provided when calling self.build_oligo()
        """
        output = []
        for subseq_obj in self.subseqs:
            if subseq_obj.is_blank() and not subseq_obj.is_optional():
                if show:
                    subseq_obj.pretty_print("id", "description")
                output.append(subseq_obj)
        return output

    def replace_seq(self, subseq_id, subseq_seq) -> None:

        """
        Replaces existing sequence of sub-sequence of given ID ('subseq_id' in
        'A_0_0' or ('A', 0, 0) or SODId formats) with given sequence ('subseq_seq').
        No return value.
        """
        subseq_obj = self.get_subseq(subseq_id)
        if subseq_obj:
            subseq_obj.sequence = subseq_seq

    def add_subseq(self, sod_subseq, pole, distance=None, index_start=1):
        """
        Adds new SODSubSeq object to sub-sequences. If distance isn't specified, the sub-sequence
        will be tagged on the end of the specified pole.
        :param sod_subseq: SODSubSeq object
        :param pole: int, describing whether this sequence should be on the 3' (3) or 5' (5) end
        :param distance: int, if provided, will shift all sub-sequences at and further from
        specified distance further from main by one position, and replace the sub-sequence
        originally at the specified distance with 'sod_subseq'
        :param index_start: int, number from which to count distance
        :return: None
        """
        subseqs = sort_sodsubseqs(self.filter_id(pole=pole))
        new_distance = distance if distance else len(subseqs)

        # update 'distance' attribute of every sub-sequence after 'distance' (var)
        if distance:
            for subseq_obj in subseqs[distance-index_start:]:
                subseq_obj.id().distance += 1

        # create new SODId and assign everything
        sod_subseq.id = SODId((self.id(), pole, new_distance))
        sod_subseq.id.set_oligo(self)
        sod_subseq.id.set_subseq(sod_subseq)
        self.subseqs.append(sod_subseq)

    def combine_oligos(self, sod_oligo) -> None:
        """
        Accepts an SODOligo object WITH A MAIN SEQUENCE (i.e. must have x_0_0 sub-sequence) and
        elongates the reference SODOligo on both 5' and 3' ends with the 5' and 3' components
        of the input SODOligo object
        :param sod_oligo: SODOligo object.
        :return: None
        """
        three_prime = sod_oligo.get_3prime_design()
        five_prime = sod_oligo.get_5prime_design()
        for subseq_obj in sort_sodsubseqs(three_prime):
            self.add_subseq(subseq_obj, 3)
        for subseq_obj in sort_sodsubseqs(five_prime):
            self.add_subseq(subseq_obj, 5)

    def build_oligo(self, *seqs, sss=False) -> str or SSSSeq:

        """
        Accepts sequences (as many as each design calls for, as (<seq_id>, <sequence>)
        tuples) and returns an oligo sequence built from them using the structure given
        by self.get_design(). seq_id should be provided in string format (e.g. 'A_0_0').

        E.g. self.build_oligos(('A_0_0', 'aaa'),
                               ('A_3_3', 'ttt')) --> 'BiotinLINKERSEQUENCEaaaLINKERSEQUENCEROLLINGCIRCLEHYBRIDISATIONSEQ1tttROLLINGCIRCLEHYBRIDISATIONSEQ2'

        :param seqs: iterable of (<seq_id>, <sequence>) tuple pairs. Format:
                self.build_oligos(('A_0_0', 'aaa'), ('A_3_3', 'ttt'))
            Optionally, if source is known for each sub-sequence, 'seqs' may be provided in
            (<seq_id>, <seq_sequence>, <path_to_seq_source>) format.
        :param sss: bool. Default = False. Set to True to return SSSSeq object.
        :return: str (if sss=False), SSSSeq object (if sss=True). Output is an oligo sequence
        built from them using the structure given by self.get_design(). seq_id should be
        provided in string format (e.g. 'A_0_0').
        """
        # corral input sequences into dictionary for easy lookup, index by tupified
        #   version of subseq id
        inputs = {}
        for seq in seqs:
            if isinstance(seq, Seq):
                seq_id, seq_seq, seq_src = seq.id, seq.seq, seq.src
            else:
                seq_id, seq_seq = seq[0:2]
                seq_src = None if len(seq) < 3 else seq[2]
            inputs[str(seq_id)] = {"id": seq_id,
                                   "seq": seq_seq,
                                   "src": seq_src}

        # iterate through all components of oligonucleotide and substitute with input
        #   sequences where necessary/possible
        output_main = ''
        output_subseqs = []
        curr_pos = order = 0
        oligo_design = self.get_design()

        # for each SODSubSeq object in oligo_design, create SSSSubSeq obj and add sequence to 'output_main'
        for subseq_obj in oligo_design:

            # get the three main fields of the current SODSubSeq (subseq_obj): id, sequence, and source
            curr_id = subseq_obj.id()
            curr_sequence: str = subseq_obj.sequence()
            curr_source = None if (not subseq_obj.oligo() or not subseq_obj.oligo().file()) else \
                subseq_obj.oligo().file().fname()

            # if current sub-sequence is a blank sub-sequence
            if subseq_obj.is_blank():

                # access the corresponding input sequence and replace current sequence fields
                #    with corresponding values in fields of input (id, sequence, source)
                try:
                    curr_id = inputs[str(curr_id)]["id"]
                    curr_sequence = inputs[str(curr_id)]["seq"]
                    curr_source = inputs[str(curr_id)]["src"]
                    
                # if none of the input sequences are for the current blank field
                except KeyError:
                    # only raise error if blank sub-sequence is NOT OPTIONAL, break out of function
                    if not subseq_obj.is_optional():
                        raise Exception("No sequence found for compulsory sub-sequence id '{}', description '{}'.".format(str(subseq_obj.id()),
                                                                                                                          subseq_obj.description()))
                    # set current sub-sequence to '' if optional and no sequence provided
                    curr_sequence = ''
                    curr_source = None

            # add the current sub-sequence to the final output sequence and create SSSSubSeq object
            output_main += curr_sequence
            output_subseqs.append(SSSSubSeq(curr_sequence, order, curr_pos,
                                            subseq_obj.description(),
                                            str(curr_id),
                                            src_path=curr_source))

            # increment relevant counters
            curr_pos += len(curr_sequence)
            order += 1
            
        return output_main if not sss else SSSSeq(output_main, ssssubseqs=output_subseqs)


class SODSubSeq:

    def __init__(self, subseq_id, sequence, description):

        self.sodid: SODId = subseq_id
        
        self.subseq_sequence: str = sequence
        self.subseq_description: str = description

        # SODOligo object to which this subseq belongs. Set with self.set_oligo()
        self.sodoligo = None

    def id(self):
        """
        :return: SODId object
        """
        return self.sodid

    def oligo_id(self):
        return self.id().oligo_id

    def pole(self):
        return self.id().pole

    def distance(self):
        return self.id().distance

    def oligo(self) -> SODOligo:
        return self.sodoligo

    def sequence(self) -> str:
        return self.subseq_sequence

    def description(self) -> str:
        return self.subseq_description

    def set_oligo(self, parent_oligo) -> None:
        self.sodoligo = parent_oligo

    # blanks are indicated by "<>" (compulsory) or "<<>>" (optional) in sequence field
    #   these subsequences are not part of the template and need to be provided when
    #   generating oligo sequences
    # TODO: allow default sequences if none provided (e.g. "<AAAA>" will output "AAAA" if
    #   no sequence was provided during oligo generation.)
    def is_blank(self) -> bool:
        return self.sequence()[0] == '<' and self.sequence()[-1] == '>'

    # optional blanks are indicated by "<<>>" in sequence field
    #   these subsequences are not required by template and may be optionally provided when
    #   generating oligo sequences
    def is_optional(self) -> bool:
        return self.sequence()[:2] == "<<" and self.sequence()[-2:] == ">>"

    def flatten(self, *fields) -> str:
        if not fields:
            fields = ('i', 's', 'd')
        to_flatten = []
        to_flatten += [self.id()] if contains_any(fields, "ID", "id", 'i') else []
        to_flatten += [self.sequence()] if contains_any(fields, "Sequence", "sequence", "seq", 's') else []
        to_flatten += [self.description()] if contains_any(fields, "Description", "description", "desc", 'd') else []
        return join_ele('\t', *to_flatten)
    
    def pretty_print(self, *fields) -> None:
        print(self.flatten(*fields))


class SODId:

    # receives ID in string format (e.g. 'A_0_0') or tuple format (e.g. ('A', 0, 0))
    def __init__(self, struct_id):

        # converts struct_id to tuple if string input is given
        if isinstance(struct_id, str):
            struct_id = re.findall(re.compile("(.+?)_(\d+?)_(\d+?)"), struct_id)[0]

        # unpack struct_id into oligo_id, pole, and distance
        self.oligo_id, self.pole, self.distance = struct_id

        # convert relevant values to integers
        self.pole = int(self.pole)
        self.distance = int(self.distance)

        # SODOligo object to which this subseq belongs. Set with self.set_oligo()
        self.sodoligo = None
        # SODSubSeq object to which this subseq belongs. Set with self.set_subseq()
        self.subseq = None

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.string()

    def __eq__(self, other):
        return self is other or\
               (isinstance(other, str) and other == self.string()) or\
               (isinstance(other, tuple) and other == self.tup())

    def string(self) -> str:
        return join_ele('_', *(self.tup()))

    def tup(self) -> tuple:
        return self.oligo_id, self.pole, self.distance

    def set_oligo(self, parent_oligo) -> None:
        self.sodoligo = parent_oligo

    def set_subseq(self, subseq) -> None:
        self.subseq = subseq

    def oligo(self) -> SODOligo:
        return self.sodoligo


######################
##       MISC       ##
##   SORTING FUNC   ##
######################

def sort_sodsubseqs(*sodsubseqs, reverse=False) -> list:
    """
    Sorts any number of SODSubSeq objects
    :param sodsubseqs: SODSubSeq objects. Accepts any quantity
    :param reverse: bool. Whether to reverse sorting order
    :return: list. List of sorted SODSubSeqs
    """
    return list(sorted(sodsubseqs, key=lambda subseq_obj: subseq_obj.id().tup(), reverse=reverse))
