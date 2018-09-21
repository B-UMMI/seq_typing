#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
get_stx_db.py - Gets STX sequences from virulencefinder_db to produce a STX subtyping DB
<https://github.com/B-UMMI/seq_typing/modules/>

Copyright (C) 2018 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: August 22, 2018

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import sys
import argparse
import os
import time

from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

from itertools import product

import utils as utils


version = '1.0'


def extend_ambiguous_dna(seq):
    """
    Return list of all possible sequences given an ambiguous DNA input
    From https://stackoverflow.com/questions/27551921/how-to-extend-ambiguous-dna-sequence

    Parameters
    ----------
    seq : str
        String with DNA sequence containing the ambiguity IUPAC codes

    Returns
    -------
    all_possible_sequences : list
        List with all different possible sequences. If memory does not allow to get all possible sequences, return None.
    """

    d = Seq.IUPAC.IUPACData.ambiguous_dna_values

    try:
        all_possible_sequences = list(map(''.join, product(*map(d.get, seq))))
    except MemoryError:
        print('ERROR: Memory Error'
              '\n')
        all_possible_sequences = None

    return all_possible_sequences


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 get_stx_db.py"')

    parser = argparse.ArgumentParser(prog='get_stx_db.py',
                                     description='Gets STX sequences from virulencefinder_db to produce a STX subtyping'
                                                 ' DB',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the sequences will be stored (default: ./)',
                                         required=False, default='.')

    args = parser.parse_args()

    start_time = time.time()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Get virulencefinder_db
    url = 'https://bitbucket.org/genomicepidemiology/virulencefinder_db.git'
    virulencefinder_db = os.path.join(args.outdir, 'virulence_db', '')
    run_successfully, _, _ = utils.runCommandPopenCommunicate(['git', 'clone', url, virulencefinder_db],
                                                              False, None, True)
    _, commit, _ = utils.runCommandPopenCommunicate(['git', '-C', virulencefinder_db, 'log',
                                                     '--pretty=format:"%h"', '-n', '1'], True, 15, True)

    # Get STX sequences
    stx_seq = {}
    # stx_seq_write = []
    allowed_chars = set(Seq.IUPAC.IUPACData.unambiguous_dna_letters)
    with open(os.path.join(args.outdir,
                           'virulence_db.virulence_ecoli.commit_{commit}.problematic_sequences.tab'.format(
                               commit=commit)), 'wt', newline='\n') as writer:
        for seq in SeqIO.parse(os.path.join(virulencefinder_db, 'virulence_ecoli.fsa'), 'fasta'):
            if seq.id.lower().startswith('stx'):
                subtype = seq.id.split(':')
                if len(subtype) == 4:
                    if seq.id[:4] not in stx_seq:
                        stx_seq[seq.id[:4]] = []

                    subtype = subtype[0][:4] + subtype[3]  # Define subtype
                    # if subtype not in stx_seq[seq_name[3]]:
                    #     stx_seq[seq_name[3]][subtype] = []
                    seq.description = ''  # To avoid description to be print in outfile

                    # For sequences with IUPAC codes, use one possible sequence based on the one with the codes
                    if not set(seq.seq.upper()).issubset(allowed_chars):
                        # print(seq.id, set(seq.seq.upper()))
                        all_possible_sequences = extend_ambiguous_dna(seq.seq.upper())
                        if all_possible_sequences is not None:
                            seq = SeqRecord(Seq.Seq(all_possible_sequences[0], generic_dna),
                                            id='{seq_name}:IUPAC_codes_removed'.format(seq_name=seq.id),
                                            description='')  # Change the sequence
                        else:
                            writer.write('\t'.join([seq.id, 'Memory Error (too much IUPAC codes)']))
                            continue

                    seq.id = '{seq_name}:seqTyping_{subtype}'.format(seq_name=seq.id, subtype=subtype)
                    stx_seq[seq.id[:4]].append(seq)
                    # stx_seq_write.append(seq)

    # Write files
    for gene, seqs in stx_seq.items():
        with open(os.path.join(args.outdir,
                               'virulence_db.virulence_ecoli.commit_{commit}.{gene}_subtyping.seq_typing.fasta'.format(
                                   commit=commit, gene=gene)), 'wt', newline='\n') as writer:
            _ = SeqIO.write(seqs, writer, "fasta")

    # print(len(stx_seq))
    # for gene, subtype_dict in stx_seq.items():
    #     print(gene, len(subtype_dict))
    #     for subtype, seqs in subtype_dict.items():
    #         print(subtype, len(seqs))

    _ = utils.runTime(start_time)


if __name__ == "__main__":
    main()
