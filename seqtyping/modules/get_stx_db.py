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
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import IUPACData

from itertools import product

try:
    import modules.utils as utils
except ImportError:
    try:
        from seqtyping.modules import utils as utils
    except ImportError:
        import utils


version = '2.0.0'


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

    d = IUPACData.ambiguous_dna_values

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
                                         help='Path to the directory where the sequences will be stored',
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
    stx_seq_seq = {}
    # stx_seq_write = []
    allowed_chars = set(IUPACData.unambiguous_dna_letters)
    with open(os.path.join(args.outdir,
                           'virulence_db.virulence_ecoli.commit_{commit}.problematic_sequences.tab'.format(
                               commit=commit)), 'wt', newline='\n') as writer:
        for seq in SeqIO.parse(os.path.join(virulencefinder_db, 'virulence_ecoli.fsa'), 'fasta'):
            if seq.id.lower().startswith('stx'):
                stx_dimer = seq.id[:4].lower()

                if len(seq.seq) > 800:
                    stx_subunit = stx_dimer + "A"
                elif 400 > len(seq.seq) > 100:
                    stx_subunit = stx_dimer + "B"
                else:
                    print(f"Undetermined stx_subunit for {seq.id}")
                    continue

                fields = seq.id.split(':')
                if len(fields[0]) >= 5:
                    if stx_subunit not in stx_seq:
                        stx_seq[stx_subunit] = []
                    if stx_subunit not in stx_seq_seq:
                        stx_seq_seq[stx_subunit] = []

                    _subtype = fields[0][4]
                    if _subtype in ["-", " ", "_", ":"]:
                        continue

                    '''
                    Jani
                    
                    After spending what seemed to be an endless amount of hours trying to solve the STEC stx subtype
                    mystery I've come to the following conclusion. For the platform we need to combine in the target db
                    stx2a,  stx2c and  stx2d as one subtype called stx2acd. This is due to the fact that all of these
                    subtypes are the most potent ones to cause HUS and cannot be separated from each other by the
                    methods in use right now.
                    '''


                    subtype = _subtype
                    if stx_dimer == 'stx2' and _subtype in ['a', 'c', 'd']:
                        subtype = 'acd'
                    else:
                        subtype = _subtype

                    subtype_str = stx_subunit + subtype  # Define subtype
                    # if subtype not in stx_seq[seq_name[3]]:
                    #     stx_seq[seq_name[3]][subtype] = []
                    seq.description = ''  # To avoid description to be print in outfile

                    # For sequences with IUPAC codes, use one possible sequence based on the one with the codes
                    seq_seq = seq.seq.upper()
                    if not set(seq_seq).issubset(allowed_chars):
                        # print(seq.id, set(seq.seq.upper()))
                        all_possible_sequences = extend_ambiguous_dna(seq_seq)
                        if all_possible_sequences is not None:
                            seq = SeqRecord(Seq.Seq(all_possible_sequences[0]),
                                            id='{seq_name}:IUPAC_codes_removed'.format(seq_name=seq.id),
                                            description='')  # Change the sequence
                        else:
                            writer.write('\t'.join([seq.id, 'Memory Error (too much IUPAC codes)']))
                            continue

                    if seq_seq not in stx_seq_seq[stx_subunit] or not any(seq_seq in seq_db for seq_db in stx_seq_seq[stx_subunit]):
                        seq_id = seq.id.replace(" ", "_")
                        seq.id = f"{seq_id}:seqTyping_{subtype_str}"
                        stx_seq[stx_subunit].append(seq)
                        # stx_seq_write.append(seq)
                        stx_seq_seq[stx_subunit].append(seq_seq)

    # Write files
    gene_counter = 1
    for gene, seqs in sorted(stx_seq.items()):
        with open(os.path.join(
                args.outdir,
                f"virulence_db.virulence_ecoli.commit_{commit}.{gene_counter}_{gene}_subtyping.seq_typing.fasta"
            ), "wt", newline="\n") as writer:
            _ = SeqIO.write(seqs, writer, "fasta")
        gene_counter += 1

    # print(len(stx_seq))
    # for gene, subtype_dict in stx_seq.items():
    #     print(gene, len(subtype_dict))
    #     for subtype, seqs in subtype_dict.items():
    #         print(subtype, len(seqs))

    _ = utils.run_time('get_stx_db.py', start_time)


if __name__ == "__main__":
    main()
