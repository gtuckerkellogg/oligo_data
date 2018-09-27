# import sys
# import os

# from oligo_errors import *
# from oligo_gen import *
from oligo_sod import *
from oligo_sss import *
from oligo_seq import *


# Version 3
# created 2018_09_26


# Accepts IDs of subseqs to fill ('subseqs_id'), sequences to fill with ('sequences'),
#   and a SODOligo object ('design_oligo'). Generates oligos based on 'design_oligo' and
#   writes generated oligos to .sss file.
# INPUTS
#   design: SODOligo object specifying oligo design
#   subseqs_id: list/tuple of sub-sequence IDs to fill (E.g. ('A_0_0', 'A_3_2'))
#   sequences: list/tuple of sub-sequence sequences. Two-level.
#       - Length of inner level should match length of 'subseqs_id'
#       - User must ensure that order of seqeunces matches order of subseq ID in 'subseqs'
#       - E.g. [('MAIN1_for_A_0_0', 'BARCODE1_for_A_3_2'),
#               ('MAIN2_for_A_0_0', 'BARCODE2_for_A_3_2'),
#               ('MAIN3_for_A_0_0', 'BARCODE3_for_A_3_2')]
#   out_fname: name of .sss to be written
#   out_dir: path to directory to which output .sss file is to be written
# OUTPUTS
#   Writes .sss file. Returns nothing.
def generate_oligos(design_oligo: SODOligo, subseqs_id, sequences,
                    out_fname: str, out_dir: str=''):

    compulsory_subseqs = design_oligo.compulsory_subseqs()

    # for SODSubSeq objects of compulsory sequences
    for compulsory_subseq in compulsory_subseqs:
        # if ID of compulsory subseq not found in user-provided list of subseq IDs
        if not str(compulsory_subseq.id()) in subseqs_id:
            raise Exception("Compulsory subsequence '{}' needs to be filled.".format(str(compulsory_subseq.id())))

    # create ((<subseq id>, <sequence>), ..., (<subseq id>, <sequence>)) tuple for each
    #   oligonucleotide to be generated
    mapped_seqs = tuple(map(lambda seqs: tuple(zip(subseqs_id, *transpose(seqs))),
                            sequences))

    output_sss = SSSFile(out_fname)
    for subseqs in mapped_seqs:
        output_sss.add_seq(design_oligo.build_oligo(*subseqs, sss=True))

    output_sss.write(out_fname, out_dir=out_dir)

    return


class OligoBuilder:

    def __init__(self, sodoligo):

        self.sodoligo = sodoligo
        self.main_subseq = sodoligo.get_0()[0]

        self.blank_sodsubseqs = {}  # dictionary of SeqFile obj indexed by str(SODSubSeq ID)
        for sodsubseq in self.sodoligo.blank_subseqs():
            self.blank_sodsubseqs[str(sodsubseq.id())] = None

    def oligo(self) -> SODOligo:
        return self.sodoligo

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

    def compulsory_seqs_filled(self):
        unfilled_subseqs = []
        for subseq in self.oligo().compulsory_subseqs():
            if self.blank_sodsubseqs[str(subseq.id())] is None:
                unfilled_subseqs.append(str(subseq.id()))
        if not unfilled_subseqs:
            return True
        print("The following compulsory sub-sequences were not provided: {}".format(join_ele(', ',
                                                                                             unfilled_subseqs)))
        return False

    # INPUTS
    #   name (str): name of SSSFile object to be created. Will be used during writing.
    #   out_dir (str): path to directory to write file. Optional if 'write' is set to False
    #   write (bool): indicates whether to write newly created SSSFile object into .sss file
    def create_oligos(self, name, out_dir='', write=True) -> SSSFile:

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
        # based on mainseq_len if x_0_x exists in list of blank subseqs, else number of sequences
        # in shortest file provided
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
