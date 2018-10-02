import re
import yaml

from oligo_gen import *
from oligo_sss import *
from oligo_seq import *


##################
##    READER    ##
##   FUNCTION   ##
##################

def read_sodyaml(fname: str):
    """
    Accepts a .yaml file for oligo designs that contains information on subsequences
    of relevant oligonucleotides for a given experiment type/amplification strategy
    and reads it into a SODFile object.
    """
    f = yaml.load(open(fname, 'r'))

    # check file type
    metadata = f["metadata"]
    if metadata.get("file type", None) != "SOD":
        while True:
            usr_in = input("'SOD' not indicated as 'file type'. Proceed anyway? Y/N\t")
            if usr_in == 'N' or usr_in == 'n':
                print("'read_sodyaml' aborted for file '{}'.".format(fname))
                return
            elif usr_in == 'Y' or usr_in == 'y':
                break

    # get data of oligos
    oligos_data = {}
    for oligo, poles_data in f["oligos"].items():
        oligos_data[oligo] = poles_data

    return SODFile(fname, oligos_data, metadata)


#####################
##   DESIGN FILE   ##
##     CLASSES     ##
#####################

class SODFile:

    def __init__(self, fname: str, oligos_data: dict, metadata: dict):

        self.name = fname
        self.metadata = metadata
        self.oligos = {}
        for oligo, poles_data in oligos_data.items():
            self.oligos[oligo] = SODOligo(oligo_id=oligo,
                                          data_yaml=poles_data,
                                          sodfile=self)

    ''' GETTER FUNCTIONS (BASIC) '''

    def fname(self) -> str:
        """
        :return: file name
        """
        return self.name

    def get_oligo(self, oligo_id):
        """
        :param oligo_id: SODOligo ID (str)
        :return: SODOligo object
        """
        return self.oligos[oligo_id]

    def get_all_subseqs(self):
        subseqs = list(map(lambda oligo_obj: oligo_obj.subseqs, self.oligos.values()))
        output = []
        for subseqs_lst in subseqs:
            output.extend(subseqs_lst)
        return output

    def get_subseq(self, subseq_id_alias):
        """
        :param subseq_id_alias: SODSubSeq ID, as string ('A_0_0') or tuple (('A', 0, 0)) or SODId object or alias (str)
        :return: SODSubSeq object that has the given ID
        """
        for subseq_obj in self.get_all_subseqs():
            if subseq_obj.is_id_alias(subseq_id_alias):
                return subseq_obj
        raise Exception("Subsequence '{}' not found in '{}'.".format(subseq_id_alias, self.fname()))

    ''' GETTER FUNCTIONS (PROCESSED) '''

    def oligo_ids(self) -> list:
        """
        :return: lst of str. List of oligo IDs.
        """
        return list(self.oligos.keys())

    def num_oligos(self) -> int:
        return len(self.oligos)

    ''' BOOLEAN FUNCTIONS '''

    def valid_oligo_id(self, oligo_id) -> bool:
        return oligo_id in self.oligo_ids()

    ''' MODIFIER FUNCTIONS '''

    def replace_sequence(self, subseq_id_alias, subseq_seq) -> None:
        subseq_obj = self.get_subseq(subseq_id_alias)
        subseq_obj.replace_sequence(subseq_seq)

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def yamlfy(self) -> dict:
        """
        Converts data in SODFile to dictionaries and lists according to YAML format
        :return: dictionary of oligos + their data
        """
        # call SODOligo.yamlfy method as well to build final YAML dict
        yamlfied = {"metadata": self.metadata}
        for oligo_id, oligo_obj in self.oligos.items():
            yamlfied[oligo_id] = oligo_obj.yamlfy()
        return yamlfied

    def write(self, fname, out_dir='') -> None:
        """
        Writes SODFile object to YAML file.
        """

        make_dir(dir_path=out_dir, check=False)
        f = open(os.path.join(out_dir, fname), "w+")
        data_yaml = self.yamlfy()

        yaml.dump(data_yaml, f, default_flow_style=False)

        f.close()

        return


class SODOligo:

    def __init__(self, oligo_id=None, data_yaml=None, sodfile=None, description=None):

        self.sodfile = sodfile
        self.oligo_id = oligo_id
        self.subseqs = []
        self.description = description

        if data_yaml:
            for pole, poles_data in data_yaml.items():
                if pole == "description":
                    self.description = poles_data
                    continue
                for distance, subseq_data in poles_data.items():
                    self.subseqs.append(SODSubSeq(sodoligo=self,
                                                  pole=pole,
                                                  distance=distance,
                                                  data_yaml=subseq_data))

    ''' GETTER FUNCTIONS (BASIC) '''

    def id(self) -> None or str:
        return self.oligo_id

    def file(self) -> SODFile:
        return self.sodfile

    ''' GETTER FUNCTIONS (PROCESSED) '''

    def filter_id(self, pole=None, distance=None) -> list:

        """
        :return: list. List of SODSubSeq objects of subsequences at a given pole
        (5' (pole = 5) or 3' (pole = 3) or none (pole = 0)) and/or distance from 0.
        Both 'pole' and 'distance' variables only support integer values.
        Not sorted in any order.
        """
        return [sodsubseq for sodsubseq in self.subseqs if
                ((pole is None or sodsubseq.pole() == pole) and
                 (distance is None or sodsubseq.distance() == distance))]

    def get_subseq(self, subseq_id_alias):

        """
        :param subseq_id_alias: subseq_id in 'A_0_0' or ('A', 0, 0) or SODId formats.
        :return: SODSubSeq object if found, else None
        """
        for subseq_obj in self.subseqs:
            if subseq_obj.is_id_alias(subseq_id_alias):
                return subseq_obj
        print("{} not found in oligo {}.".format(str(subseq_id_alias), self.id()))

    def get_3prime_oligos(self) -> list:
        """
        :return: list. List of SODSubSeq objects of subsequences at the 3' end.
        Not sorted in any order.
        """
        return self.filter_id(pole=3)

    def get_5prime_oligos(self) -> list:

        """
        :return: list. List of SODSubSeq objects of subsequences at the 5' end.
        Not sorted in any order.
        """
        return self.filter_id(pole=5)

    def get_0(self) -> list:

        """
        :return: list. List of SODSubSeq object of subsequence at 0.
        """
        return self.filter_id(pole=0)

    def get_design(self) -> list:

        """
        Returns a list of SODSubSeq objects of subsequences in 5' --> 3' order.
        """
        five_prime = sort_sodsubseqs(*self.get_5prime_oligos(), reverse=True)
        three_prime = sort_sodsubseqs(*self.get_3prime_oligos())
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

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def set_sodfile(self, sodfile) -> None:
        """
        Set SODFile object to which this SODOligo belongs
        :param sodfile: SODFile
        :return: None
        """
        self.sodfile = sodfile

    def replace_sequence(self, subseq_id_alias, subseq_seq) -> None:

        """
        Replaces existing sequence of sub-sequence of given ID ('subseq_id' in
        'A_0_0' or ('A', 0, 0) or SODId formats) with given sequence ('subseq_seq').
        No return value.
        """
        subseq_obj = self.get_subseq(subseq_id_alias)
        if subseq_obj:
            subseq_obj.replace_sequence(subseq_seq)

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
            for subseq_obj in subseqs[distance - index_start:]:
                subseq_obj.id().distance += 1

        # create new SODId and assign everything
        sod_subseq.set_id(SODId((self.id(), pole, new_distance)))
        sod_subseq.id().set_oligo(self)
        sod_subseq.id().set_subseq(sod_subseq)
        self.subseqs.append(sod_subseq)

    def combine_oligos(self, sod_oligo) -> None:
        """
        Accepts an SODOligo object WITH A MAIN SEQUENCE (i.e. must have x_0_0 sub-sequence) and
        elongates the reference SODOligo on both 5' and 3' ends with the 5' and 3' components
        of the input SODOligo object
        :param sod_oligo: SODOligo object.
        :return: None
        """
        three_prime = sod_oligo.get_3prime_oligos()
        five_prime = sod_oligo.get_5prime_oligos()
        for subseq_obj in sort_sodsubseqs(three_prime):
            self.add_subseq(subseq_obj, 3)
        for subseq_obj in sort_sodsubseqs(five_prime):
            self.add_subseq(subseq_obj, 5)

    ''' DISPLAY FUNCITONS '''

    def pretty_print_subseqs(self, *fields) -> None:
        """
        Prints the values of subsequences at specified fields  in pretty-ish format
        :param fields: str. Fields to print. Options: "id" ('i'), "sequence" ("seq", 's'), "description" ("desc", 'd')
        :return: None
        """

        header = []
        header += ["ID"] if contains_any(fields, "ID", "id", 'i') else []
        header += ["Sequence"] if contains_any(fields, "Sequence", "sequence", "seq", 's') else []
        header += ["Description"] if contains_any(fields, "Description", "description", "desc", 'd') else []
        print(join_ele('\t', *header))

        for subseq_obj in self.subseqs:
            subseq_obj.pretty_print(*header)

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def yamlfy(self) -> dict:
        oligos = {0: self.get_0(),
                  3: self.get_3prime_oligos(),
                  5: self.get_5prime_oligos()}

        yamlfied = {}
        for pole, subseqs in oligos.items():
            if subseqs:
                yamlfied[pole] = {}
                for subseq_obj in subseqs:
                    yamlfied[pole][subseq_obj.distance()] = subseq_obj.yamlfy()

        return yamlfied

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
            curr_source = None if not self.file() else self.file().fname()

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
                        raise Exception(
                            "No sequence found for compulsory sub-sequence id '{}', description '{}'.".format(
                                subseq_obj.id(),
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

    def __init__(self, oligo=None, pole=None, distance=None, alias=None, data_yaml=None, sodoligo=None):

        self.sodid: SODId = SODId(oligo_id=(None if (not oligo and not sodoligo) else
                                            sodoligo.id()),
                                  pole=pole, distance=distance, alias=data_yaml.get("alias", alias),
                                  sodsubseq=self)

        self.subseq_sequence = data_yaml.get("sequence", None)
        self.subseq_description = data_yaml.get("description", None)
        self.subseq_identical = data_yaml.get("identical", None)
        self.reverse_complement = data_yaml.get("reverse complement", None)

        # SODOligo object to which this subseq belongs. Set with self.set_oligo()
        self.sodoligo = sodoligo

    ''' GETTER FUNCTIONS (BASIC) '''

    def id(self):
        """
        :return: SODId object
        """
        return self.sodid

    def oligo_id(self) -> str:
        return self.id().oligo_id

    def oligo(self) -> SODOligo:
        return self.sodoligo

    def pole(self) -> int:
        return self.id().pole

    def distance(self) -> int:
        return self.id().distance

    def alias(self) -> str:
        return self.id().alias()

    def description(self) -> str:
        return self.subseq_description

    ''' GETTER FUNCTIONS (PROCESSED) '''

    def sequence(self) -> str:
        """
        Returns self.subseq_sequence directly if it is not None.
        Else, dynamically looks up sequence if a sub-sequence ID or a sequence is found in self.complement.
        Else, raises Exception.
        :return: str. ATGC sequence of self.
        """
        if self.subseq_sequence and self.valid_sequence():
            return self.subseq_sequence
        elif self.reverse_complement:
            ref = self.oligo().file().get_subseq(self.reverse_complement)
            if not ref or not ref.valid_sequence():
                raise Exception("No sequence complement reference found for {}.".format(self.id()))
            return ref.sequence_rcomp()
        elif self.subseq_identical:
            ref = self.oligo().file().get_subseq(self.subseq_identical)
            if not ref or not ref.valid_sequence():
                raise Exception("No identical sequence reference found for {}.".format(self.id()))
            return ref.sequence()
        else:
            raise Exception("No valid sequence or complement found for {}.".format(self.id()))

    def sequence_rcomp(self) -> str:
        """
        Generates reverse complement of self's sequence.
        :return: str. Reverse complement of self.
        """
        lookup_table = {'A': 'T', 'a': 't',
                        'T': 'A', 't': 'a',
                        'C': 'G', 'c': 'g',
                        'G': 'C', 'g': 'c'}
        return join_ele('', *map(lambda char: lookup_table[char], self.sequence()[::-1]))

    ''' BOOLEAN FUNCTIONS '''

    # blanks are indicated by "<>" (compulsory) or "<<>>" (optional) in sequence field
    #   these subsequences are not part of the template and need to be provided when
    #   generating oligo sequences
    # TODO: allow default sequences if none provided (e.g. "<AAAA>" will output "AAAA" if
    #   no sequence was provided during oligo generation.)
    def is_blank(self) -> bool:
        return False if not self.subseq_sequence else\
            (self.subseq_sequence[0] == '<' and self.subseq_sequence[-1] == '>')

    # optional blanks are indicated by "<<>>" in sequence field
    #   these subsequences are not required by template and may be optionally provided when
    #   generating oligo sequences
    def is_optional(self) -> bool:
        return self.is_blank and\
               self.subseq_sequence[:2] == "<<" and self.subseq_sequence[-2:] == ">>"

    def is_id_alias(self, id_alias) -> bool:
        return self.id() == id_alias or self.alias() == id_alias

    def valid_sequence(self) -> bool:
        valid_char = ['a', 't', 'g', 'c', 'A', 'T', 'G', 'C']
        return False if not self.subseq_sequence else\
            True if (self.is_blank()) else\
            (len(tuple(filter(lambda char: char in valid_char, self.subseq_sequence))) ==
             len(self.subseq_sequence))

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def set_id(self, new_sodid) -> None:
        self.sodid = new_sodid

    def set_oligo(self, parent_oligo) -> None:
        """
        Assigns SODOligo object to self. Updates self's SODId object to reflect change.
        :param parent_oligo: SODOligo object.
        :return: None
        """
        self.sodoligo = parent_oligo
        self.sodid.set_oligo(parent_oligo)

    def set_rcomplement(self, subseq_id):
        self.reverse_complement = subseq_id

    def replace_sequence(self, sequence):
        self.subseq_sequence = sequence

    ''' DISPLAY FUNCTIONS '''

    def pretty_print(self, *fields) -> None:
        print(self.flatten(*fields))

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def flatten(self, *fields) -> str:
        if not fields:
            fields = ('i', 's', 'd')
        to_flatten = []
        to_flatten += [self.id()] if contains_any(fields, "ID", "id", 'i') else []
        to_flatten += [self.sequence()] if contains_any(fields, "Sequence", "sequence", "seq", 's') else []
        to_flatten += [self.description()] if contains_any(fields, "Description", "description", "desc", 'd') else []
        return join_ele('\t', *to_flatten)

    def yamlfy(self) -> dict:
        """
        Converts data into dictionary for YMAL file format.
        {"sequence": <sequence info>,
         "reverse complement": <reverse complement reference>,
         "description": <description>}
        Empty fields are removed.
        :return: dict.
        """
        yamlfied = {"alias": self.id().alias(),
                    "sequence": self.subseq_sequence,
                    "reverse complement": self.reverse_complement,
                    "description": self.description()}
        to_del = []
        for field, value in yamlfied.items():
            if value is None:
                to_del.append(field)
        for field in to_del:
            del yamlfied[field]
        return yamlfied


class SODId:

    # receives ID broken down into oligo_id, pole, and distance, or
    # in string format (e.g. 'A_0_0') or tuple format (e.g. ('A', 0, 0))
    def __init__(self, *alt_format, oligo_id=None, pole=None, distance=None, alias=None, sodsubseq=None, sodoligo=None):

        self.oligo_id = oligo_id
        self.pole = pole
        self.distance = distance
        self.subseq_alias = alias

        # converts struct_id to tuple if string input is given
        if alt_format:
            alt_format = alt_format[0]
            # if is alias
            if isinstance(alt_format, str) and alt_format[0] == '_':
                self.subseq_alias = alt_format
            # else if is an alternative form of ID (e.g. ('A', 0, 0) or 'A_0_0')
            else:
                if isinstance(alt_format, str):
                    alt_format = re.findall(re.compile("(.+?)_(\d+?)_(\d+?)"), alt_format)[0]
                self.oligo_id, self.pole, self.distance = alt_format

        # convert relevant values to integers
        self.pole = int(self.pole)
        self.distance = int(self.distance)

        # SODSubSeq object to which this subseq belongs. Set with self.set_subseq()
        self.subseq = sodsubseq
        # SODOligo object to which this subseq belongs. Set with self.set_oligo()
        self.sodoligo = sodoligo

    ''' OVER-RIDING FUNCTIONS '''

    def __repr__(self):
        return str(self)

    def __str__(self):
        return self.string()

    def __eq__(self, other):
        return self is other or \
               (isinstance(other, str) and other == self.string()) or \
               (isinstance(other, tuple) and other == self.tup())

    ''' GETTER FUNCTIONS (BASIC) '''

    def oligo(self) -> SODOligo:
        return self.sodoligo

    def alias(self) -> str:
        return self.subseq_alias

    ''' CONVERSION FUNCTIONS '''

    def string(self) -> str:
        return join_ele('_', *(self.tup()))

    def tup(self) -> tuple:
        return self.oligo_id, self.pole, self.distance

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def set_oligo(self, parent_oligo) -> None:
        self.sodoligo = parent_oligo
        self.oligo_id = parent_oligo.id()

    def set_subseq(self, subseq) -> None:
        self.subseq = subseq

    def set_alias(self, new_alias) -> None:
        self.subseq_alias = new_alias


######################
##       MISC       ##
##   SORTING FUNC   ##
######################

def sort_sodsubseqs(*sodsubseqs, reverse=False) -> list:
    """
    Sorts any number of SODSubSeq objects by SODId (in order of oligo, pole, distance)
    :param sodsubseqs: SODSubSeq objects. Accepts any quantity
    :param reverse: bool. Whether to reverse sorting order
    :return: list. List of sorted SODSubSeqs
    """
    return list(sorted(sodsubseqs, key=lambda subseq_obj: subseq_obj.id().tup(), reverse=reverse))
