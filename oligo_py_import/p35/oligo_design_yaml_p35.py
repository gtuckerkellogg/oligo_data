from .oligo_sodyaml_p35 import *
from .oligo_sss_p35 import *
from ..oligo_seq import *


# Version 3
# created 2018_09_26


class OligoBuilder:

    def __init__(self, sodoligo):

        self.sodoligo = sodoligo
        self.main_subseq = sodoligo.get_0()[0]

        self.blank_sodsubseqs = {}  # dictionary of SeqFile obj indexed by str(SODSubSeq ID)
        for sodsubseq in self.sodoligo.blank_subseqs():
            self.blank_sodsubseqs[str(sodsubseq.id())] = None

    ''' GETTER FUNCTIONS (BASIC) '''

    def oligo(self) -> SODOligo:
        return self.sodoligo

    def get_seqfile(self, subseq_id):
        """
        :param subseq_id: accepts an SODSubSeq ID (string, tuple, or SODId object are okay)
        :return: SeqFile object that was assigned to the given sub-sequence, or None if no
            object was assigned
        """
        for blank_subseq in self.oligo().blank_subseqs():
            curr_subseq_id = blank_subseq.id()
            if curr_subseq_id == subseq_id:
                return self.blank_sodsubseqs[str(curr_subseq_id)]
        print("No SeqFile object was found for sub-sequence ID {}.".format(str(subseq_id)))

    ''' BOOLEAN FUNCTIONS '''

    def compulsory_seqs_filled(self) -> bool:
        unfilled_subseqs = []
        for subseq in self.oligo().compulsory_subseqs():
            if self.blank_sodsubseqs[str(subseq.id())] is None:
                unfilled_subseqs.append(str(subseq.id()))
        if not unfilled_subseqs:
            return True
        print("The following compulsory sub-sequences were not provided: {}".format(join_ele(', ',
                                                                                             unfilled_subseqs)))
        return False

    ''' SETTER AND MODIFIER FUNCTIONS '''

    def set_main(self, subseq_id) -> None:
        self.main_subseq = self.oligo().get_subseq(subseq_id)

    def add_seqfile(self, subseq_id, seqfile) -> None:
        """
        If the subseq_id already has a SeqFile object assigned to it, it will be overwritten
        with the newly given SeqFile object ('seqfile')
        :param subseq_id: string of ID of a SODSubSeq (e.g. 'A_0_0') to assign this file to
        :param seqfile: SeqFile object pointing to file with sequences
        :return: None
        """
        if subseq_id not in self.blank_sodsubseqs:
            raise Exception("{} is not a valid sub-sequence ID for Oligo {} from file {}.".format(subseq_id,
                                                                                                  self.oligo().id(),
                                                                                                  self.oligo().file().fname()))
        self.blank_sodsubseqs[subseq_id] = seqfile
        return

    ''' DATA REPACKAGING AND GENERATION FUNCTIONS '''

    def create_oligos(self, name: str, out_dir: str ='', write: bool =True) -> SSSFile:
        """
        Generates oligos.
        :param name: str. Name of SSSFile to generate. Will be used during writing.
        :param out_dir: str. Directory to write SSSFile object to if write is set to True. Default: <working directory>
        :param write: bool. SSSFile object to be written to file if set to True. Default: True
        :return: SSSFile object.
        """

        # checks if compulsory sub-sequences are filled.
        if not self.compulsory_seqs_filled():
            raise Exception("Unable to generate oligos because compulsory sub-sequences were not provided.")

        sssfile = SSSFile(name)

        # lists of SODSubSeq IDs and SeqFile object, where indices for corresponding
        # SODSubSeq ID - SeqFile object pairs are the same
        blank_ids = [subseq.id() for subseq in self.oligo().blank_subseqs()
                     if self.blank_sodsubseqs[str(subseq.id())] is not None]
        blank_seqfiles = [self.blank_sodsubseqs[str(subseq_id)]
                          for subseq_id in blank_ids]

        # get number of sequences provided for main oligo subsequence (ID: x_0_x). If it isn't
        # a blank subseq, set this to None
        mainseq_len = None if not self.main_subseq\
            else self.get_seqfile(self.main_subseq.id()).num_seqs()
        # set number of times to loop through blank_subseqs_objs to generate new sequences
        #   based on mainseq_len if x_0_x exists in list of blank subseqs, else number of
        #   sequences in shortest file provided
        range_ref = (mainseq_len if mainseq_len
                     else min(map(lambda seqfile_obj: seqfile_obj.num_seqs(), blank_seqfiles))) - 1

        for i in range(range_ref):

            # corral sequences into format acceptable by SODOligo.build_oligo()
            zipped_seqs = tuple(zip(blank_ids,
                                    list(map(lambda seqfile: seqfile.seq_at(i).seq,
                                             blank_seqfiles)),
                                    list(map(lambda seqfile: seqfile.seq_at(i).src,
                                             blank_seqfiles))))

            # generate new oligo and add newly created SSSSeq obj to sssfile
            sssfile.add_seq(self.oligo().build_oligo(*zipped_seqs, sss=True))

        # write sequences.
        if write:
            sssfile.write(name, out_dir=out_dir)

        return sssfile
