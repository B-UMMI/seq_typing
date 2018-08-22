#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import sys
from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

from itertools import product
from collections import Counter


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

    counter_all = Counter(seq)
    counter_iupac = 0
    counter_iupac_scored = 0
    for letter, value in counter_all.items():
        if letter not in list(Seq.IUPAC.IUPACData.unambiguous_dna_letters):
            counter_iupac += value
            counter_iupac_scored += value * len(d[letter])

    all_possible_sequences = None
    if counter_iupac_scored / float(len(seq)) < 0.002:
        try:
            all_possible_sequences = list(map(''.join, product(*map(d.get, seq))))
        except MemoryError:
            print('ERROR: Memory Error\n'
                  'Frequency of IUPAC codes: {counter_iupac}\n'
                  'Frequency of IUPAC codes scored: {counter_iupac_scored}\n'
                  '\n'.format(counter_iupac=counter_iupac / float(len(seq)),
                              counter_iupac_scored=counter_iupac_scored / float(len(seq))))
            all_possible_sequences = None

    return all_possible_sequences


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling using "python3 check_sequences.py"')

    fasta_file = 'GenotypesDENV_14-05-18.fasta'

    # Get sequences
    sequences = {}
    allowed_chars = set('ATGC')
    with_gap = 0
    for seq in SeqIO.parse(fasta_file, 'fasta'):
        if not seq.id.lower().split('|')[4].startswith('subtype:'):
            if 'no_subtype' not in sequences:
                sequences['no_subtype'] = []
            print('NO SUBTYPE', seq.id)
            sequences['no_subtype'].append(seq)
        else:
            if not set(seq.seq.upper()).issubset(set('-')):
                seq = SeqRecord(Seq.Seq(str(seq.seq).replace('-', ''), generic_dna),
                                id='{seq_name}|gaps_removed'.format(seq_name=seq.id),
                                description='')  # Change the sequence
                with_gap += 1

            if not set(seq.seq.upper()).issubset(allowed_chars):
                if 'iupac_code' not in sequences:
                    sequences['iupac_code'] = []
                # print('IUPAC CODES', seq.id, set(seq.seq.upper()))
                all_possible_sequences = extend_ambiguous_dna(seq.seq.upper())
                if all_possible_sequences is not None:
                    seq = SeqRecord(Seq.Seq(all_possible_sequences[0], generic_dna),
                                    id='{seq_name}|IUPAC_codes_removed'.format(seq_name=seq.id),
                                    description='')  # Change the sequence
                    sequences['iupac_code'].append(seq)
                else:
                    if 'impossible_iupac_code' not in sequences:
                        sequences['impossible_iupac_code'] = []
                    sequences['impossible_iupac_code'].append(seq)
            else:
                seq.id = '{seq_name}|seqTyping_{subtype}'.format(seq_name=seq.id,
                                                                 subtype=seq.id.split('|')[4].lstrip('Subtype:'))
                seq.description = ''  # To avoid description to be print in outfile
                if len(seq.id.split('|')[4].split('_')) > 1:
                    if 'with_underscore' not in sequences:
                        sequences['with_underscore'] = []
                    # print('WITH UNDERSCORE', seq.id)
                    sequences['with_underscore'].append(seq)
                else:
                    if 'normal' not in sequences:
                        sequences['normal'] = []
                    sequences['normal'].append(seq)

    print('gap', with_gap)
    for type_seq, seqs in sequences.items():
        print(type_seq, len(seqs))

    # Write files
    with open(fasta_file + '.no_gap_iupac.fasta', 'wt', newline='\n') as writer:
        _ = SeqIO.write(sequences['normal'] + sequences['iupac_code'], writer, 'fasta')

    # Get similar sequences
    # Exact matches
    sequences_exact = {}
    counter = 0
    for seq in SeqIO.parse(fasta_file + '.no_gap_iupac.fasta', 'fasta'):
        if seq.seq in sequences_exact:
            sequences_exact[seq.seq].append(seq.id)
        else:
            sequences_exact[seq.seq] = [seq.id]
        counter += 1
    print('# Sequences to analyse: {n}'.format(n=counter))
    print('# Different sequences: {n}'.format(n=len(sequences_exact)))

    # Contained sequences
    sequences_contained = {}
    counter = 0
    for seq_exact, ids_exact in sequences_exact.items():
        not_in_seq_contained = True
        for seq_contained, ids_contained in sequences_contained.items():
            if seq_exact in seq_contained or seq_contained in seq_exact:
                if len(seq_exact) <= len(seq_contained):
                    # First ids_contained to be coherent with the sequence used as key (seq_contained)
                    sequences_contained[seq_contained] = ids_contained + ids_exact
                else:
                    del sequences_contained[seq_contained]
                    # First ids_exact to be coherent with the sequence used as key (seq_exact)
                    sequences_contained[seq_exact] = ids_exact + ids_contained
                not_in_seq_contained = False
                counter += 1
                break
        if not_in_seq_contained:
            sequences_contained[seq_exact] = ids_exact
    print('# Sequences contained: {n}'.format(n=counter))
    print('# Unique sequences (excluding contained sequences): {n}'.format(n=len(sequences_contained)))

    # Write unique sequences
    sequences_to_write = []
    for seq in SeqIO.parse(fasta_file + '.no_gap_iupac.fasta', 'fasta'):
        # Can use 0 in sequences_contained[record.seq][0] due to coherence established above
        if seq.seq in sequences_contained and seq.id == sequences_contained[seq.seq][0]:
            sequences_to_write.append(seq)
    with open(fasta_file + '.no_gap_iupac.fasta.unique_sequences.fasta', 'wt', newline='\n') as writer:
        _ = SeqIO.write(sequences_to_write, writer, 'fasta')

    # Link between unique sequences and other sequences
    with open(fasta_file + '.no_gap_iupac.fasta.unique_sequences.relation.tab', 'wt') as writer:
        writer.write('\t'.join(['#seq_used', 'similar_sequences']) + '\n')
        for seq in SeqIO.parse(fasta_file + '.no_gap_iupac.fasta.unique_sequences.fasta', 'fasta'):
            writer.write('\t'.join([sequences_contained[seq.seq][0]] +
                                   [';'.join(sequences_contained[seq.seq][1:]) if len(sequences_contained[seq.seq]) > 1
                                    else ''])
                         + '\n')


if __name__ == "__main__":
    main()


'''
gap 4026
normal 3581
iupac_code 364
impossible_iupac_code 81
no_subtype 7

# Sequences to analyse: 3945
# Different sequences: 3819
# Sequences contained: 1
# Unique sequences (excluding contained sequences): 3818
'''
