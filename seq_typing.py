#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
seq_typing.py - Determine which reference sequence is more likely to be present
in a given sample
<https://github.com/B-UMMI/seq_typing/>

Copyright (C) 2017 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: June 26, 2017

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

import argparse
import os
import time
import sys

import modules.utils as utils
import modules.run_rematch as run_rematch
import modules.parse_results as parse_results

version = '0.1'


def parse_config(config_file):
    config = {'length_extra_seq': None, 'minimum_depth_presence': None, 'minimum_depth_call': None,
              'minimum_gene_coverage': None}

    with open(config_file, 'rtU') as reader:
        field = None
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                line = line.split(' ')[0]
                if line.startswith('#'):
                    line = line[1:].split(' ')[0]
                    field = line
                else:
                    if field is not None:
                        if field in ['length_extra_seq', 'minimum_depth_presence', 'minimum_depth_call',
                                     'minimum_gene_coverage']:
                            line = int(line)
                            if field == 'minimum_gene_coverage':
                                if line < 0 or line > 100:
                                    sys.exit('minimum_gene_coverage in config file must be an integer between 0 and'
                                             ' 100')
                        config[field] = line
                        field = None

    for field in config:
        if config[field] is None:
            sys.exit('{} in config file is missing'.format(field))

    return config


def get_fasta_config(species):
    """
    Get the reference fasta file and config file for the species provided

    Parameters
    ----------
    species : list
        List with strings that correspond to the species name, e.g. ('escherichia', 'coli')

    Returns
    -------
    fasta : list
        Sorted list of fasta files to be used as references
    config : str
        File path to config file
    """

    fasta = []
    config = None

    file_path = os.path.abspath(__file__)
    species_folder = os.path.join(os.path.dirname(file_path), 'serotyping_reference_sequences', '_'.join(species), '')

    files = [f for f in os.listdir(species_folder) if not f.startswith('.') and
             os.path.isfile(os.path.join(species_folder, f))]
    for file_found in files:
        file_found = os.path.join(species_folder, file_found)
        if file_found.endswith('.fasta'):
            fasta.append(file_found)
        elif file_found.endswith('.config'):
            config = str(file_found)

    return sorted(fasta), config


def get_species_allowed():
    """
    Get the name of the species with reference sequences provided for serotyping

    Returns
    -------
    species : list
        List with species names, e.g. ['escherichia coli', 'streptococcus agalactiae']
    """

    file_path = os.path.abspath(__file__)
    serotyping_folder = os.path.join(os.path.dirname(file_path), 'serotyping_reference_sequences', '')
    species = [d.replace('_', ' ') for d in os.listdir(serotyping_folder) if not d.startswith('.') and
               os.path.isdir(os.path.join(serotyping_folder, d))]
    return species


def include_rematch_dependencies_path():
    command = ['which', 'rematch.py']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
    if run_successfully:
        rematch_script = stdout.splitlines()[0]
        utils.setPATHvariable(False, rematch_script)
        return rematch_script
    else:
        sys.exit('ReMatCh not found in the PATH')


def clean_header(header, problematic_characters):
    new_header = header
    if any(x in header for x in problematic_characters):
            for x in problematic_characters:
                new_header = new_header.replace(x, '_')
    return header, new_header


def parse_reference(reference, problematic_characters):
    reference_dict = {}
    headers_correspondence = {}
    with open(reference, 'rtU') as reader:
        header = None
        sequence = ''
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if line.startswith('>'):
                    if header is not None:
                        reference_dict[header] = sequence
                    original_header, new_header = clean_header(line[1:], problematic_characters)
                    if new_header in headers_correspondence:
                        sys.exit('Possible conflicting sequence header in {reference} file:'.format(reference=reference) + '\n' +
                                 '{original_header} header might be the same as {first_header} header after problematic characters ({problematic_characters}) replacement (new header: {new_header})'.format(original_header=original_header, first_header=headers_correspondence[new_header], problematic_characters=problematic_characters, new_header=new_header))
                    header = str(new_header)
                    headers_correspondence[header] = str(original_header)
                    sequence = ''
                else:
                    sequence += line.replace(' ', '').upper()
        if len(sequence) > 0:
            reference_dict[header] = sequence
    return reference_dict, headers_correspondence


def rename_duplicated_headers(references_headers, reference, reference_dict, headers_correspondence,
                              problematic_characters):
    renamed_reference_dict, renamed_headers_correspondence, headers_changed = {}, {}, []
    for ref in references_headers:
        if any(x in references_headers[ref].keys() for x in reference_dict):
            for header in reference_dict:
                if header in references_headers[ref]:
                    original_header, new_header = clean_header('_'.join([header,
                                                                         os.path.basename(reference),
                                                                         'SeqTyping']),
                                                               problematic_characters)
                    renamed_reference_dict[new_header] = reference_dict[header]
                    renamed_headers_correspondence[new_header] = headers_correspondence[header]
                    headers_changed.append(header)
                    utils.Bcolors_print('{header} sequence header from {reference} file was found already'
                                        ' in {ref} file. Therefore it was renamed'
                                        ' to {new_header}'.format(header=header, reference=os.path.basename(reference),
                                                                  ref=os.path.basename(ref),
                                                                  new_header=new_header), 'WARNING')
    for header in reference_dict:
        if header not in headers_changed:
            renamed_reference_dict[header] = reference_dict[header]
            renamed_headers_correspondence[header] = headers_correspondence[header]
    return renamed_reference_dict, renamed_headers_correspondence


def write_sequence(reference_dict, writer):
    for header, sequence in reference_dict.items():
        writer.write('>' + header + '\n')
        writer.write('\n'.join(utils.chunkstring(sequence, 80)) + '\n')


def prepare_references(references, map_ref_together, references_dir):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    references_headers = {}
    references_files = []

    if map_ref_together:
        ref_file = os.path.join(references_dir, 'references_together.fasta')
        writer = open(ref_file, 'wt')
        references_files.append(ref_file)

    for reference in references:
        reference_dict, headers_correspondence = parse_reference(reference, problematic_characters)
        reference_dict, headers_correspondence = rename_duplicated_headers(references_headers, reference,
                                                                           reference_dict, headers_correspondence,
                                                                           problematic_characters)

        references_headers[reference] = dict(headers_correspondence)

        if not map_ref_together:
            ref_file = os.path.join(references_dir, os.path.basename(reference))
            writer = open(ref_file, 'wt')
            references_files.append(ref_file)

        write_sequence(reference_dict, writer)
        writer.flush()

        if not map_ref_together:
            writer.close()

    if map_ref_together:
        writer.close()

    return references_files, references_headers


def main():
    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 seq_typing.py"')

    parser = argparse.ArgumentParser(prog='seq_typing.py',
                                     description='Determine which reference sequence is more likely to be present in a'
                                                 ' given sample',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version', version=str('%(prog)s v' + version))

    parser_required = parser.add_argument_group('Required options')
    parser_required.add_argument('-f', '--fastq', nargs='+', action=utils.required_length((1, 2), '--fastq'), type=argparse.FileType('r'), metavar=('/path/to/input/file.fq.gz'), help='Path to single OR paired-end fastq files. If two files are passed, they will be assumed as being the paired fastq files', required=True)

    parser_reference = parser.add_mutually_exclusive_group()
    parser_reference.add_argument('-r', '--reference', nargs='+', type=argparse.FileType('r'),
                                  metavar='/path/to/reference_sequence.fasta',
                                  help='Fasta file containing reference sequences. If more than one file is passed, a'
                                       ' reference sequence for each file will be determined. Give the files name in'
                                       ' the same order that the type must be determined.')
    parser_reference.add_argument('-s', '--species', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                  help='Name of the species with reference sequences provided together with %(prog)s'
                                       ' for serotyping',
                                  action=utils.arguments_choices_words(get_species_allowed(), '--species'))

    parser_optional_general = parser.add_argument_group('General facultative options')
    parser_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/', help='Path to the directory where the information will be stored', required=False, default='.')
    parser_optional_general.add_argument('-j', '--threads', type=int, metavar='N', help='Number of threads to use', required=False, default=1)
    parser_optional_general.add_argument('--mapRefTogether', action='store_true', help='Map the reads against all references together')
    parser_optional_general.add_argument('--typeSeparator', type=str, metavar='_', help='Last single character separating the general sequence header from the last part containing the type', required=False, default='_')
    parser_optional_general.add_argument('--extraSeq', type=int, metavar='N', help='Sequence length added to both ends of target sequences (usefull to improve reads mapping to the target one) that will be trimmed in ReMatCh outputs', required=False, default=0)
    parser_optional_general.add_argument('--minCovPresence', type=int, metavar='N', help='Reference position minimum coverage depth to consider the position to be present in the sample', required=False, default=5)
    parser_optional_general.add_argument('--minCovCall', type=int, metavar='N', help='Reference position minimum coverage depth to perform a base call', required=False, default=10)
    parser_optional_general.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
                                         help=argparse.SUPPRESS, required=False, default=0.6)
    # parser_optional_general.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
    #                                      help='Minimum relative frequency of the dominant allele coverage depth'
    #                                           ' (value between [0, 1]). Positions with lower values will be'
    #                                           ' considered as having multiple alleles (and will be coded as N)',
    #                                      required=False, default=0.6)
    parser_optional_general.add_argument('--minGeneCoverage', type=int, metavar='N',
                                         help='Minimum percentage of target reference sequence covered to consider a'
                                              ' sequence to be present (value between [0, 100])',
                                         required=False, default=60)
    parser_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N', help=argparse.SUPPRESS,
                                         required=False, default=80)
    # parser_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N',
    #                                      help='Minimum percentage of identity of reference sequence covered to'
    #                                           ' consider a gene to be present (value between [0, 100]). One INDEL'
    #                                           ' will be considered as one difference',
    #                                      required=False, default=80)
    parser_optional_general.add_argument('--doNotRemoveConsensus', action='store_true', help='Do not remove ReMatCh consensus sequences')
    parser_optional_general.add_argument('--debug', action='store_true', help='DeBug Mode: do not remove temporary files')
    parser_optional_general.add_argument('--beginning', action='store_true', help='Start %(prog)s from the beggining')
    parser_optional_general.add_argument('--notClean', action='store_true', help='Do not remove intermediate files')

    args = parser.parse_args()

    if args.minGeneCoverage is not None and (args.minGeneCoverage < 0 or args.minGeneCoverage > 100):
        parser.error('--minGeneCoverage should be a value between [0, 100]')
    # if args.minGeneIdentity is not None and (args.minGeneIdentity < 0 or args.minGeneIdentity > 100):
    #     parser.error('--minGeneIdentity should be a value between [0, 100]')
    if args.reference is None and args.species is None:
        parser.error('--reference or --species must be provided')

    start_time = time.time()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Start logger
    logfile, time_str = utils.start_logger(args.outdir)

    script_path = utils.general_information(logfile, version, args.outdir, time_str)
    del script_path
    print('\n')

    rematch_script = include_rematch_dependencies_path()

    args.fastq = [fastq.name for fastq in args.fastq]

    if args.reference is not None:
        args.reference = [reference.name for reference in args.reference]
    else:
        args.reference, config = get_fasta_config(args.species)
        config = parse_config(config)
        args.extraSeq = config['length_extra_seq']
        args.minCovPresence = config['minimum_depth_presence']
        args.minCovCall = config['minimum_depth_call']
        args.minGeneCoverage = config['minimum_gene_coverage']

        print('\n'
              'Settings that will be used:\n'
              '    reference: {reference}\n'
              '    extraSeq: {extraSeq}\n'
              '    minCovPresence: {minCovPresence}\n'
              '    minCovCall: {minCovCall}\n'
              '    minGeneCoverage: {minGeneCoverage}\n'
              '\n'.format(reference=args.reference, extraSeq=args.extraSeq, minCovPresence=args.minCovPresence,
                          minCovCall=args.minCovCall, minGeneCoverage=args.minGeneCoverage))

    folders_2_remove = []

    # Create modules pickles folder
    pickles_folder = os.path.join(args.outdir, 'pickles', '')
    if not os.path.isdir(pickles_folder):
        os.makedirs(pickles_folder)
    folders_2_remove.append(pickles_folder)

    # Prepare references
    references_dir = os.path.join(args.outdir, 'references', '')
    if not os.path.isdir(references_dir):
        os.makedirs(references_dir)
    folders_2_remove.append(references_dir)
    references_files, references_headers = prepare_references(args.reference, args.mapRefTogether, references_dir)

    # Run ReMatCh
    pickleFile = os.path.join(pickles_folder, 'rematch_module.pkl')
    if not os.path.isfile(pickleFile) or args.beginning:
        runtime, references_results, module_dir = run_rematch.run_rematch(rematch_script, args.outdir, references_files, args.fastq, args.threads, args.extraSeq, args.minCovPresence, args.minCovCall, args.minFrequencyDominantAllele, args.minGeneCoverage, args.minGeneIdentity, args.debug, args.doNotRemoveConsensus)
        folders_2_remove.append(module_dir)
        utils.saveVariableToPickle([references_results, module_dir], pickleFile)
    else:
        print('ReMatCh module already run')
        references_results, module_dir = utils.extractVariableFromPickle(pickleFile)
        folders_2_remove.append(module_dir)

    # Parse ReMatCh results
    pickleFile = os.path.join(pickles_folder, 'parse_results.pkl')
    if not os.path.isfile(pickleFile) or args.beginning:
        seq_type, seq_type_info, probable_results, improbable_results = parse_results.parse_results(references_results, args.reference, references_headers, args.outdir, args.minGeneCoverage, args.typeSeparator)
        utils.saveVariableToPickle([seq_type, seq_type_info, probable_results, improbable_results], pickleFile)
    else:
        print('Results parser module already run')
        seq_type, seq_type_info, probable_results, improbable_results = utils.extractVariableFromPickle(pickleFile)
        parse_results.write_reports(args.outdir, seq_type, seq_type_info, probable_results, improbable_results)

    if not args.notClean and not args.debug:
        for folder in folders_2_remove:
            utils.removeDirectory(folder)

    time_taken = utils.runTime(start_time)
    del(time_taken)


if __name__ == "__main__":
    main()
