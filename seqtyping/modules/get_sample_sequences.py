#!/usr/bin/env python3

# -*- coding: utf-8 -*-

import os, sys, argparse, time


try:
    from __init__ import __version__ as _VERSION_
except ImportError:
    from seqtyping.__init__ import __version__ as _VERSION_


_PROGRAM_NAME_ = os.path.basename(__file__)


def get_sample_sequences(seqtyping_outdir, seqtyping_refdir, outdir='.'):

    seqtyping_report_types_tsv = [f for f in os.listdir(seqtyping_outdir) if 
            f.startswith('seq_typing.') and f.endswith('.report_types.tab') and
            os.path.isfile(os.path.join(seqtyping_outdir, f))]

    if len(seqtyping_report_types_tsv) != 1:
        raise RuntimeError(
                'It was not possible to find only one seq_typing.*.report_types.tab. {c} were found: {l}'.format(
                c=len(seqtyping_report_types_tsv), l=seqtyping_report_types_tsv))
    else:
        seqtyping_report_types_tsv = os.path.join(seqtyping_outdir, seqtyping_report_types_tsv[0])

    selected_sequences = {}
    with open(seqtyping_report_types_tsv, 'rt') as reader:
        for line in reader:
            if not line.startswith('#'):
                line = line.strip().split('\t')
                if line[0] == 'selected':  # #sequence_type
                    selected_sequences[os.path.basename(line[1])] = {  # reference_file
                                'type': line[2],
                                'sequence_id': line[3],
                                'sequence': None
                            }


    seqtyping_new_allele_dir = [d for d in os.listdir(seqtyping_outdir) if 
            d == 'new_allele' and
            os.path.isdir(os.path.join(seqtyping_outdir, d))]

    if len(seqtyping_new_allele_dir) > 1:
        raise RuntimeError(
                'It was not possible to find only one new_allele directory. {c} were found: {l}'.format(
                c=len(seqtyping_new_allele_dir), l=seqtyping_new_allele_dir))
    elif len(seqtyping_new_allele_dir) == 1:
        seqtyping_new_allele_dir = os.path.join(seqtyping_outdir, seqtyping_new_allele_dir[0])
    else:
        seqtyping_new_allele_dir = None

    if seqtyping_new_allele_dir is not None:
        x = [d for d in os.listdir(seqtyping_new_allele_dir) if 
                os.path.isdir(os.path.join(seqtyping_new_allele_dir, d))]
        for reference_file in x:
            if reference_file in selected_sequences:
                reference_dir = os.path.join(seqtyping_new_allele_dir, reference_file)
                y = [f for f in os.listdir(reference_dir) if 
                        f == selected_sequences[reference_file]['type'] + '.fasta' and
                        os.path.isfile(os.path.join(reference_dir, f))]
                for type_fasta in y:
                    with open(os.path.join(reference_dir, type_fasta), 'rt') as reader:
                        selected_sequences[reference_file]['sequence'] = []
                        for line in reader:
                            selected_sequences[reference_file]['sequence'].append(line.strip())


    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    for reference_file in selected_sequences:
        if selected_sequences[reference_file]['sequence'] is None:
            real_reference_file = os.path.join(seqtyping_refdir, reference_file)
            sequence = []
            with open(real_reference_file, 'rt') as reader:
                found_sequence = False
                for line in reader:
                    line = line.rstrip('\r\n')
                    if line.startswith('>') and not found_sequence:
                        if line == '>' + selected_sequences[reference_file]['sequence_id']:
                            sequence.append(line)
                            found_sequence = True
                    elif line.startswith('>') and found_sequence:
                        break
                    elif not line.startswith('>') and found_sequence:
                        sequence.append(line)
            if len(sequence) > 0:
                selected_sequences[reference_file]['sequence'] = sequence
        
        with open(os.path.join(outdir, reference_file), 'wt', newline='\n') as writer:
            print('\n'.join(selected_sequences[reference_file]['sequence']), file=writer)


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 {}"'.format(_PROGRAM_NAME_))

    parser = argparse.ArgumentParser(prog=_PROGRAM_NAME_,
                                     description='Gets sample sequences (now only by reference file)',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + _VERSION_))


    parser_required_general = parser.add_argument_group('Required options')
    parser_required_general.add_argument('-s', '--seqtyping_dir', type=str, metavar='/path/to/input/seqtyping_outdir/',
                                       help='Path to seq_typing output directory',
                                       required=True)
    parser_required_general.add_argument('-r', '--reference_dir', type=str, metavar='/path/to/input/seqtyping_reference_dir/',
                                       help='Path to seq_typing reference directory',
                                       required=True)

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                         help='Path to the directory where the sequences will be stored',
                                         required=False, default='.')

    args = parser.parse_args()

    start_time = time.time()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    get_sample_sequences(seqtyping_outdir=os.path.abspath(args.seqtyping_dir),
                         seqtyping_refdir=os.path.abspath(args.reference_dir),
                         outdir=args.outdir)


if __name__ == "__main__":
    main()
