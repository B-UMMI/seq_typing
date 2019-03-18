#!/usr/bin/env python3

# -*- coding: utf-8 -*-

"""
seq_typing.py - Determine which reference sequence is more likely to be present
in a given sample
<https://github.com/B-UMMI/seq_typing/>

Copyright (C) 2019 Miguel Machado <mpmachado@medicina.ulisboa.pt>

Last modified: January 10, 2019

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
from pkg_resources import resource_filename
from shutil import copyfile as shutil_copyfile

try:
    from __init__ import __version__

    import modules.utils as utils
    import modules.run_rematch as run_rematch
    import modules.parse_results as parse_results
    import modules.run_blast as run_blast
except ImportError:
    from seqtyping.__init__ import __version__

    from seqtyping.modules import utils as utils
    from seqtyping.modules import run_rematch as run_rematch
    from seqtyping.modules import parse_results as parse_results
    from seqtyping.modules import run_blast as run_blast


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

    file_path = os.path.abspath(os.path.realpath(__file__))
    species_folder = os.path.join(os.path.dirname(file_path), 'reference_sequences', '_'.join(species), '')

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

    file_path = os.path.abspath(os.path.realpath(__file__))
    serotyping_folder = os.path.join(os.path.dirname(file_path), 'reference_sequences', '')
    species = [d.replace('_', ' ') for d in os.listdir(serotyping_folder) if not d.startswith('.') and
               os.path.isdir(os.path.join(serotyping_folder, d))]
    return species


def include_rematch_dependencies_path():
    original_rematch = None
    command = ['which', 'rematch.py']
    run_successfully, stdout, stderr = utils.runCommandPopenCommunicate(command, False, None, False)
    if run_successfully:
        original_rematch = stdout.splitlines()[0]

    resource_rematch = None
    try:
        resource_rematch = resource_filename('ReMatCh', 'rematch.py')
    except ModuleNotFoundError:
        resource_rematch = original_rematch
    else:
        print('\n'
              'Using ReMatCh "{resource_rematch}" via "{original_rematch}"\n'.format(resource_rematch=resource_rematch,
                                                                                     original_rematch=original_rematch))

    if resource_rematch is not None:
        utils.setPATHvariable(False, resource_rematch)
    else:
        sys.exit('ReMatCh not found in the PATH')

    return resource_rematch


def rename_duplicated_headers(references_headers, reference, reference_dict, headers_correspondence,
                              problematic_characters):
    renamed_reference_dict, renamed_headers_correspondence, headers_changed = {}, {}, []
    for ref in references_headers:
        if any(x in references_headers[ref].keys() for x in reference_dict):
            for header in reference_dict:
                if header in references_headers[ref]:
                    original_header, new_header = utils.clean_header('_'.join([header,
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
        reference_dict, headers_correspondence = utils.parse_reference(reference, problematic_characters)
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


def assembly_subcommand(args):
    msg = []
    # if args.blast is None and args.org is None:
    #     msg.append('--blast or --org must be provided')
    if args.blast is not None and args.type is None:
        msg.append('With --blast option you must provide the --type')
    if args.minGeneCoverage < 0 or args.minGeneCoverage > 100:
        msg.append('--minGeneCoverage should be a value between [0, 100]')
    if args.minGeneIdentity < 0 or args.minGeneIdentity > 100:
        msg.append('--minGeneIdentity should be a value between [0, 100]')

    if len(msg) > 0:
        argparse.ArgumentParser(prog='assembly subcommand options').error('\n'.join(msg))

    if args.type == 'nucl' or args.type is None:
        utils.required_programs({'blastn': ['-version', '>=', '2.6.0']})
    elif args.type == 'prot':
        utils.required_programs({'blastp': ['-version', '>=', '2.6.0']})

    args.fasta = os.path.abspath(args.fasta[0].name)
    if args.blast is not None:
        args.blast = [os.path.abspath(blast.name) for blast in args.blast]
    else:
        args.blast, config = get_fasta_config(args.org)
        config = parse_config(config)
        if args.type != 'nucl':
            print('\n'
                  'ATTENTION: Blast DB type provided was not "nucl"\n'
                  'It was changed to "nucl"'
                  '\n')
        args.type = 'nucl'
        args.minGeneCoverage = config['minimum_gene_coverage']
        args.typeSeparator = '_'

        print('\n'
              'Settings that will be used:\n'
              '    DB reference file: {reference}\n'
              '    Blast DB type: nucl\n'
              '    minGeneCoverage: {minGeneCoverage}\n'
              '    Type separator character: {typeSeparator}'
              '\n'.format(reference=args.blast, minGeneCoverage=args.minGeneCoverage, typeSeparator=args.typeSeparator))

    folders_2_remove_all = []
    references_results_all = {}
    references_headers_all = {}
    blast_files = []
    for blast in args.blast:
        _, folders_2_remove, blast_results, blast, headers_correspondence = run_blast.run_blast(blast, args.outdir,
                                                                                                args.type, args.fasta)
        folders_2_remove_all.extend(folders_2_remove)
        references_results_all[blast] = blast_results
        references_headers_all[blast] = headers_correspondence
        blast_files.append(blast)

    return folders_2_remove_all, references_results_all, blast_files, references_headers_all


def blast_subcommand(args):
    msg = []
    if args.fasta is not None and args.type is None:
        msg.append('With --fasta option you must provide the --type')
    # if args.fasta is None and args.org is None:
    #     msg.append('--fasta or --org must be provided')

    if len(msg) > 0:
        argparse.ArgumentParser(prog='blast subcommand options').error('\n'.join(msg))

    utils.required_programs({'makeblastdb': ['-version', '>=', '2.6.0']})

    if args.fasta is not None:
        args.fasta = [os.path.abspath(fasta.name) for fasta in args.fasta]
    else:
        args.fasta, _ = get_fasta_config(args.org)
        if args.type != 'nucl':
            print('\n'
                  'ATTENTION: Blast DB type provided was not "nucl"\n'
                  'It was changed to "nucl"'
                  '\n')
        args.type = 'nucl'

        print('\n'
              'Settings that will be used:\n'
              '    fasta: {reference}\n'
              '    Blast DB type: nucl\n'
              '\n'.format(reference=args.fasta))

    utils.removeDirectory(os.path.join(args.outdir, 'pickles', ''))

    error_msg = []
    for fasta in args.fasta:
        # Create DB
        blast_db = os.path.join(args.outdir, '{blast_DB}'.format(blast_DB=os.path.basename(fasta)))
        db_exists, original_file = run_blast.check_db_exists(blast_db)
        if not db_exists and not original_file:
            db_exists = run_blast.create_blast_db(fasta, blast_db, args.type)
            if db_exists:
                print('Blast DB created for {file} in {outdir}'.format(file=fasta, outdir=args.outdir))
                # sys.exit(0)
            else:
                error_msg.append('It was not possible to create Blast DB or {}'.format(fasta))
        elif db_exists and original_file:
            error_msg.append('Blast DB already found for {file} in {outdir} as {blast_db}'.format(file=fasta,
                                                                                                  outdir=args.outdir,
                                                                                                  blast_db=blast_db))
        else:
            error_msg.append('It was found only Blast DB files or the original fasta file from which the Blast DB'
                             ' should be produced ({file}). Either include the missing files or remove the ones present'
                             ' (usually the original fasta file)'.format(file=fasta))

    if len(error_msg) == 0:
        sys.exit(0)
    else:
        sys.exit('\n'.join(error_msg))


def check_reference_exist(reference):
    """
    Check if a reference file exists, both the fasta file and the Bowtie2 index

    Parameters
    ----------
    reference : str
        Path to the reference sequences file

    Returns
    -------
    fasta_file : bool
        If the fasta file exists returns True, else False
    index_files : list
        List with Bowtie2 index files found
    pickle_file : bool
        If the pickle file for the fasta file exists returns True, else False
    """

    things_found = os.listdir(os.path.dirname(reference))

    basename_reference = os.path.basename(reference)

    fasta_file = False
    index_files = []
    pickle_file = False
    for thing in things_found:
        if thing.startswith(basename_reference):
            file_path = os.path.join(os.path.dirname(reference), thing)
            file_link_exists = os.path.isfile(file_path) or os.path.islink(file_path)
            if thing == basename_reference and file_link_exists:
                fasta_file = True
            elif thing.endswith('.bt2') and file_link_exists:
                index_files.append(file_path)
            elif thing.endswith('.pkl') and file_link_exists:
                pickle_file = True

    return fasta_file, index_files, pickle_file


def write_seq_from_sequence_dict(sequence_dict, out_fasta_file):
    """
    Write a fasta file with sequences found in sequence_dict dictionary

    Parameters
    ----------
    sequence_dict : dict
        Dictionary as follow
        sequence_dict[counter] = {'header': without_greater_sign, 'sequence': seq, 'length': seq_length_int}
    out_fasta_file : str
        Path to output fasta file

    Returns
    -------

    """
    with open(out_fasta_file, 'wt', newline='\n') as writer:
        for _, seq_info in sequence_dict.items():
            writer.write('>{header}\n'
                         '{seq}\n'.format(header=seq_info['header'],
                                          seq='\n'.join(utils.chunkstring(seq_info['sequence'], 80))))


def reads_subcommand(args):
    msg = []
    # if args.reference is None and args.org is None:
    #     argparse.ArgumentParser.error('--reference or --org must be provided')
    if args.minGeneCoverage < 0 or args.minGeneCoverage > 100:
        msg.append('--minGeneCoverage should be a value between [0, 100]')
    if args.minGeneIdentity < 0 or args.minGeneIdentity > 100:
        msg.append('--minGeneIdentity should be a value between [0, 100]')

    if len(msg) > 0:
        argparse.ArgumentParser(prog='assembly subcommand options').error('\n'.join(msg))

    rematch_script = include_rematch_dependencies_path()

    utils.required_programs({'rematch.py': ['--version', '>=', '4.0.1']})

    args.fastq = [os.path.abspath(fastq.name) for fastq in args.fastq]

    folders_2_remove = []

    references_dir = os.path.join(args.outdir, 'references', '')
    if not os.path.isdir(references_dir):
        os.makedirs(references_dir)
    folders_2_remove.append(references_dir)

    references_headers = {}

    clean_run_rematch = True

    if args.reference is not None:

        clean_run_rematch = False

        args.reference = [os.path.abspath(reference) for reference in args.reference]

        reference_files = {}

        for reference in args.reference:
            fasta_file, index_files, pickle_file = check_reference_exist(reference)
            if not fasta_file and len(index_files) == 0 and not pickle_file:
                sys.exit('Missing reference fasta file, Bowtie2 index of pickle file for {}'.format(reference))
            reference_files[reference] = fasta_file, index_files, pickle_file

        references_to_use = []

        for reference in args.reference:

            reference_files_found = reference_files[reference]

            reference_file = reference

            # Create symlink to pickle file
            if reference_files_found[2]:
                symlink = os.path.join(references_dir, os.path.basename(reference_file + '.pkl'))
                if os.path.islink(symlink):
                    os.unlink(symlink)
                os.symlink(reference_file + '.pkl', symlink)
                header_gene_list, _ = utils.extractVariableFromPickle(reference_file + '.pkl')
            else:
                new_reference_file, header_gene_list, seq_reference_dict = \
                    run_rematch.clean_headers_reference_file(reference_file=reference_file, outdir=references_dir)
                if new_reference_file != reference_file:
                    utils.Bcolors_print('WARNING: Sequences headers were renamed for {}'.format(reference), 'WARNING')
                    reference_file = new_reference_file

                pickle_file = os.path.join(references_dir, os.path.basename(reference_file) + '.pkl')
                utils.saveVariableToPickle((header_gene_list, seq_reference_dict), pickle_file)

            if reference_file == reference:

                # Create symlinks to reference file
                if reference_files_found[0]:
                    symlink = os.path.join(references_dir, os.path.basename(reference_file))
                    if os.path.islink(symlink):
                        os.unlink(symlink)
                    os.symlink(reference_file, symlink)
                    reference_file = symlink
                else:
                    # Create reference file if does not exist

                    # From index files
                    if len(reference_files_found[1]) > 0:
                        run_successfully, reference_file = \
                            run_rematch.run_bowtie_inspect(index_without_sufix=reference_file, outdir=references_dir)
                        if not run_successfully:
                            sys.exit('Something went wrong while creating the reference fasta file from Bowtie2'
                                     ' index {}'.format(reference_file))

                    # From pickle file
                    elif reference_files_found[2]:
                        _, seq_reference_dict = utils.extractVariableFromPickle(reference_file + '.pkl')
                        reference_file = os.path.join(references_dir, os.path.basename(reference_file))
                        write_seq_from_sequence_dict(seq_reference_dict, reference_file)

                # Create symlinks to index files
                if len(reference_files_found[1]) > 0:
                    for index_file in reference_files_found[1]:
                        symlink = os.path.join(references_dir, os.path.basename(index_file))
                        if os.path.islink(symlink):
                            os.unlink(symlink)
                        os.symlink(index_file, symlink)

            references_to_use.append(reference_file)
            references_headers[reference_file] = dict(header_gene_list)

        args.reference = references_to_use
    else:
        args.reference, config = get_fasta_config(args.org)

        references_to_use = []

        for reference in args.reference:
            new_reference_file, header_gene_list, seq_reference_dict = \
                run_rematch.clean_headers_reference_file(reference_file=reference, outdir=references_dir)

            if new_reference_file != reference:
                utils.Bcolors_print('WARNING: Sequences headers were renamed for {}'.format(reference), 'WARNING')
                references_to_use.append(new_reference_file)
                references_headers[new_reference_file] = dict(header_gene_list)
                pickle_file = os.path.join(references_dir, os.path.basename(new_reference_file) + '.pkl')
            else:
                symlink = os.path.join(references_dir, os.path.basename(reference))
                if os.path.islink(symlink):
                    os.unlink(symlink)
                os.symlink(reference, symlink)
                references_to_use.append(symlink)
                references_headers[symlink] = dict(header_gene_list)
                pickle_file = os.path.join(references_dir, os.path.basename(symlink) + '.pkl')

            utils.saveVariableToPickle((header_gene_list, seq_reference_dict), pickle_file)

        args.reference = references_to_use

        config = parse_config(config)
        args.extraSeq = config['length_extra_seq']
        args.minCovPresence = config['minimum_depth_presence']
        args.minCovCall = config['minimum_depth_call']
        args.minGeneCoverage = config['minimum_gene_coverage']
        args.typeSeparator = '_'

        print('\n'
              'Settings that will be used:\n'
              '    reference: {reference}\n'
              '    extraSeq: {extraSeq}\n'
              '    minCovPresence: {minCovPresence}\n'
              '    minCovCall: {minCovCall}\n'
              '    minGeneCoverage: {minGeneCoverage}\n'
              '    Type separator character: {typeSeparator}'
              '\n'.format(reference=args.reference, extraSeq=args.extraSeq, minCovPresence=args.minCovPresence,
                          minCovCall=args.minCovCall, minGeneCoverage=args.minGeneCoverage,
                          typeSeparator=args.typeSeparator))

    pickles_folder = os.path.join(args.outdir, 'pickles', '')

    # Run ReMatCh
    pickle_file = os.path.join(pickles_folder, 'rematch_module.pkl')
    if args.resume and os.path.isfile(pickle_file):
        print('ReMatCh module already run')
        references_results, module_dir = utils.extractVariableFromPickle(pickle_file)
        folders_2_remove.append(module_dir)
    else:
        _, references_results, module_dir = run_rematch.run_rematch(rematch_script, args.outdir, args.reference,
                                                                    args.fastq, args.threads, args.extraSeq,
                                                                    args.minCovPresence, args.minCovCall,
                                                                    args.minFrequencyDominantAllele,
                                                                    args.minGeneCoverage, args.minGeneIdentity,
                                                                    args.debug, args.doNotRemoveConsensus,
                                                                    args.bowtieAlgo,
                                                                    clean_run_rematch=clean_run_rematch)
        folders_2_remove.append(module_dir)
        utils.saveVariableToPickle([references_results, module_dir], pickle_file)

    return folders_2_remove, references_results, args.reference, references_headers


def index_subcommand(args):
    _ = include_rematch_dependencies_path()

    utils.required_programs({'bowtie2-build': ['--version', '>=', '2.2.9']})

    if args.reference is not None:
        args.reference = [os.path.abspath(reference.name) for reference in args.reference]
    else:
        args.reference, _ = get_fasta_config(args.org)
        print('\n'
              'Settings that will be used:\n'
              '    reference: {reference}\n'
              '\n'.format(reference=args.reference))

    for reference in args.reference:
        reference_file_outdir = os.path.join(args.outdir, os.path.basename(reference))

        new_reference_file, header_gene_list, seq_reference_dict = \
            run_rematch.clean_headers_reference_file(reference_file=reference, outdir=args.outdir)

        if new_reference_file == reference:
            # Copy fasta file to outdir
            shutil_copyfile(reference, reference_file_outdir, follow_symlinks=True)
        else:
            utils.Bcolors_print('WARNING: Sequences headers were renamed for {}'.format(reference), 'WARNING')
            # Rename files
            os.rename(new_reference_file, reference_file_outdir)

        pickle_file = os.path.join(args.outdir, reference_file_outdir + '.pkl')
        utils.saveVariableToPickle((header_gene_list, seq_reference_dict), pickle_file)

        # Create Bowtie2 index
        run_successfully = run_rematch.run_bowtie_build(reference_file=reference_file_outdir, outdir=args.outdir,
                                                        threads=args.threads)
        if not run_successfully:
            sys.exit('Something went wrong while creating Bowtie2 index for {}'.format(reference))

    exit(0)


def python_arguments(program_name, version):
    """
    Sets pythons arguments

    Parameters
    ----------
    program_name : str
        String with the name of the program
    version : str
        String with version

    Returns
    -------
    parser : argparse.ArgumentParser
        General argparse
    parser_reads : argparse.ArgumentParser.add_subparsers.add_parser
        Reads subparser
        For running the program using fastq files
    parser_index : argparse.ArgumentParser.add_subparsers.add_parser
        Index subparser
        For creating Bowtie2 index
    parser_assembly : argparse.ArgumentParser.add_subparsers.add_parser
        Assembly subparser
        For running the program using a fasta file
    parser_blast : argparse.ArgumentParser.add_subparsers.add_parser
        Blast subparser
        For creating Blast DB
    """
    parser = argparse.ArgumentParser(prog=program_name,
                                     description='Determines which reference sequence is more likely to be present in a'
                                                 ' given sample',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--version', help='Version information', action='version',
                        version='{prog} v{version}'.format(prog=parser.prog, version=version))

    subparsers = parser.add_subparsers(title='Subcommands', description='Valid subcommands', help='Additional help')

    parser_reads = subparsers.add_parser('reads', description='Run {prog} using fastq files. If running multiple'
                                                              ' samples using the same reference sequences file,'
                                                              ' consider use first "{prog} index"'
                                                              ' subcommand.'.format(prog=parser.prog),
                                         help='reads --help', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_index = subparsers.add_parser('index', description='Creates Bowtie2 index. This is useful when running the'
                                                              ' same reference sequences file for different reads'
                                                              ' dataset.',
                                         help='index --help', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_assembly = subparsers.add_parser('assembly', description='Run {prog} using a fasta file. If running'
                                                                    ' multiple samples using the same DB sequence file,'
                                                                    ' consider use first "{prog} blast"'
                                                                    ' subcommand.'.format(prog=parser.prog),
                                            help='assembly --help',
                                            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_blast = subparsers.add_parser('blast', description='Creates Blast DB. This is useful when running the same'
                                                              ' DB sequence file for different assemblies.',
                                         help='blast --help', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser_reads_required = parser_reads.add_argument_group('Required options')
    parser_reads_required.add_argument('-f', '--fastq', nargs='+', action=utils.required_length((1, 2), '--fastq'),
                                       type=argparse.FileType('r'), metavar='/path/to/input/file.fq.gz',
                                       help='Path to single OR paired-end fastq files. If two files are passed, they'
                                            ' will be assumed as being the paired fastq files',
                                       required=True)

    parser_reads_reference = parser_reads.add_mutually_exclusive_group(required=True)
    parser_reads_reference.add_argument('-r', '--reference', nargs='+', type=str,
                                        metavar='/path/to/reference/sequences.file',
                                        help='Path to reference sequences file. If Bowtie2 index was already produced,'
                                             ' only provide the file name that ends with ".1.bt2", but without this'
                                             ' termination (for example, for a Bowtie2 index'
                                             ' "/file/sequences.fasta.1.bt2", only provide "/file/sequences.fasta").'
                                             ' If no Bowtie2 index files are found, those will be created in --outdir.'
                                             ' If more than one file is passed, a type for each file will be'
                                             ' determined. Give the files name in the same order that the type must be'
                                             ' determined.')
    parser_reads_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                        help='Name of the organism with reference sequences provided together'
                                             ' with {} for typing ("seqtyping/reference_sequences/"'
                                             ' folder)'.format(parser.prog),
                                        action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_reads_optional_general = parser_reads.add_argument_group('General facultative options')
    parser_reads_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                               help='Path to the directory where the information will be stored'
                                                    ' (default: ./',
                                               required=False, default='.')
    parser_reads_optional_general.add_argument('-j', '--threads', type=int, metavar='N',
                                               help='Number of threads to use (default: 1)', required=False, default=1)
    parser_reads_optional_general.add_argument('--mapRefTogether', action='store_true',
                                               help=argparse.SUPPRESS)
    # parser_reads_optional_general.add_argument('--mapRefTogether', action='store_true',
    #                                            help='Map the reads against all references together')
    parser_reads_optional_general.add_argument('--typeSeparator', type=str, metavar='_',
                                               help='Last single character separating the general sequence header from'
                                                    ' the last part containing the type (default: _)',
                                               required=False, default='_')
    parser_reads_optional_general.add_argument('--extraSeq', type=int, metavar='N',
                                               help='Sequence length added to both ends of target sequences (usefull to'
                                                    ' improve reads mapping to the target one) that will be trimmed in'
                                                    ' ReMatCh outputs (default when not using --org: 0)',
                                               required=False, default=0)
    parser_reads_optional_general.add_argument('--minCovPresence', type=int, metavar='N',
                                               help='Reference position minimum coverage depth to consider the position'
                                                    ' to be present in the sample (default when not using --org: 5)',
                                               required=False, default=5)
    parser_reads_optional_general.add_argument('--minCovCall', type=int, metavar='N',
                                               help='Reference position minimum coverage depth to perform a base call'
                                                    ' (default when not using --org: 10)',
                                               required=False, default=10)
    parser_reads_optional_general.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
                                               help=argparse.SUPPRESS, required=False, default=0.6)
    # parser_optional_general.add_argument('--minFrequencyDominantAllele', type=float, metavar='0.6',
    #                                      help='Minimum relative frequency of the dominant allele coverage depth'
    #                                           ' (value between [0, 1]). Positions with lower values will be'
    #                                           ' considered as having multiple alleles (and will be coded as N)',
    #                                      required=False, default=0.6)
    parser_reads_optional_general.add_argument('--minGeneCoverage', type=int, metavar='N',
                                               help='Minimum percentage of target reference sequence covered to'
                                                    ' consider a sequence to be present (value between [0, 100])'
                                                    ' (default when not using --org: 60)',
                                               required=False, default=60)
    parser_reads_optional_general.add_argument('--minDepthCoverage', type=int, metavar='N',
                                               help='Minimum depth of coverage of target reference sequence to'
                                                    ' consider a sequence to be present (default: 2)',
                                               required=False, default=2)
    # parser_reads_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N', help=argparse.SUPPRESS,
    #                                            required=False, default=80)
    parser_reads_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N',
                                               help='Minimum percentage of identity of reference sequence covered to'
                                                    ' consider a gene to be present (value between [0, 100]). One INDEL'
                                                    ' will be considered as one difference (default: 80)',
                                               required=False, default=80)
    parser_reads_optional_general.add_argument('--bowtieAlgo', type=str, metavar='"--very-sensitive-local"',
                                               help='Bowtie2 alignment mode. It can be an end-to-end alignment'
                                                    ' (unclipped alignment) or local alignment (soft clipped'
                                                    ' alignment). Also, can choose between fast or sensitive'
                                                    ' alignments. Please check Bowtie2 manual for extra information:'
                                                    ' http://bowtie-bio.sourceforge.net/bowtie2/index.shtml .'
                                                    ' This option should be provided between quotes and starting with'
                                                    ' an empty space (like --bowtieAlgo " --very-fast") or using equal'
                                                    ' sign (like --bowtieAlgo="--very-fast")',
                                               required=False, default='--very-sensitive-local')
    parser_reads_optional_general.add_argument('--doNotRemoveConsensus', action='store_true',
                                               help='Do not remove ReMatCh consensus sequences')
    parser_reads_optional_general.add_argument('--debug', action='store_true',
                                               help='Debug mode: do not remove temporary files')
    parser_reads_optional_general.add_argument('--resume', action='store_true',
                                               help='Resume %(prog)s')
    parser_reads_optional_general.add_argument('--notClean', action='store_true',
                                               help='Do not remove intermediate files')

    parser_index_reference = parser_index.add_mutually_exclusive_group(required=True)
    parser_index_reference.add_argument('-r', '--reference', nargs='+', type=str,
                                        metavar='/path/to/reference/sequences.fasta',
                                        help='Path to reference sequences file. If more than one file is passed, a'
                                             ' Bowtie2 index for each file will be created.')
    parser_index_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                        help='Name of the organism with reference sequences provided'
                                             ' ("seqtyping/reference_sequences/" folder) together'
                                             ' with seq_typing.py for typing',
                                        action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_index_optional_general = parser_index.add_argument_group('General facultative options')
    parser_index_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                               help='Path to the directory where the information will be stored'
                                                    ' (default: ./)',
                                               required=False, default='.')
    parser_index_optional_general.add_argument('-j', '--threads', type=int, metavar='N', required=False,
                                               help='Number of threads to use (default: 1)', default=1)

    parser_assembly_required = parser_assembly.add_argument_group('Required options')
    parser_assembly_required.add_argument('-f', '--fasta', nargs=1, type=argparse.FileType('r'),
                                          metavar='/path/to/query/assembly_file.fasta',
                                          help='Path to fasta file containing the query sequences from which the'
                                               ' types should be assessed',
                                          required=True)

    parser_assembly_reference = parser_assembly.add_mutually_exclusive_group(required=True)
    parser_assembly_reference.add_argument('-b', '--blast', nargs='+', type=argparse.FileType('r'),
                                           metavar='/path/to/Blast/db.sequences.file',
                                           help='Path to DB sequence file. If Blast DB was already produced, only'
                                                ' provide the file that do not end with ".n*" something (do not use for'
                                                ' example /blast_db.sequences.fasta.nhr). If no Blast DB is found for'
                                                ' the DB sequence file, one will be created in --outdir. If more than'
                                                ' one Blast DB file is passed, a type for each file will be determined.'
                                                ' Give the files in the same order that the type must be determined.')
    parser_assembly_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                           help='Name of the organism with DB sequence file provided'
                                                ' ("seqtyping/reference_sequences/" folder) together'
                                                ' with seq_typing.py for typing',
                                           action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_assembly_optional_reference = parser_assembly.add_argument_group('Required option for --blast')
    parser_assembly_optional_reference.add_argument('-t', '--type', choices=['nucl', 'prot'], type=str, metavar='nucl',
                                                    help='Blast DB type (available options: %(choices)s)')

    parser_assembly_optional_general = parser_assembly.add_argument_group('General facultative options')
    parser_assembly_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                                  help='Path to the directory where the information will be stored'
                                                       ' (default: ./)',
                                                  required=False, default='.')
    parser_assembly_optional_general.add_argument('-j', '--threads', type=int, metavar='N', required=False,
                                                  help='Number of threads to use (default: 1)', default=1)
    parser_assembly_optional_general.add_argument('--typeSeparator', type=str, metavar='_',
                                                  help='Last single character separating the general sequence header'
                                                       ' from the last part containing the type (default: _)',
                                                  required=False, default='_')
    parser_assembly_optional_general.add_argument('--minGeneCoverage', type=int, metavar='N',
                                                  help='Minimum percentage of target reference sequence covered to'
                                                       ' consider a sequence to be present (value between [0, 100])'
                                                       ' (default when not using --org: 60)',
                                                  required=False, default=60)
    parser_assembly_optional_general.add_argument('--minGeneIdentity', type=int, metavar='N',
                                                  help='Minimum percentage of identity of reference sequence covered'
                                                       ' to consider a gene to be present (value between [0, 100])'
                                                       ' (default: 80)',
                                                  required=False, default=80)
    parser_assembly_optional_general.add_argument('--minDepthCoverage', type=int, metavar='N', help=argparse.SUPPRESS,
                                                  required=False, default=1)
    parser_assembly_optional_general.add_argument('--debug', action='store_true',
                                                  help='Debug mode: do not remove temporary files')
    parser_assembly_optional_general.add_argument('--resume', action='store_true', help=argparse.SUPPRESS)

    parser_blast_required = parser_blast.add_argument_group('Required options')
    parser_blast_required.add_argument('-t', '--type', choices=['nucl', 'prot'], type=str, metavar='nucl',
                                       help='Blast DB type (available options: %(choices)s)', required=True)

    parser_blast_reference = parser_blast.add_mutually_exclusive_group(required=True)
    parser_blast_reference.add_argument('-f', '--fasta', nargs='+', type=argparse.FileType('r'),
                                        metavar='/path/to/db.sequences.fasta',
                                        help='Path to DB sequence file. If more than one file is passed, a Blast DB for'
                                             ' each file will be created.')
    parser_blast_reference.add_argument('--org', nargs=2, type=str.lower, metavar=('escherichia', 'coli'),
                                        help='Name of the organism with DB sequence file provided'
                                             ' ("seqtyping/reference_sequences/" folder) together'
                                             ' with seq_typing.py for typing',
                                        action=utils.arguments_choices_words(get_species_allowed(), '--org'))

    parser_blast_optional_general = parser_blast.add_argument_group('General facultative options')
    parser_blast_optional_general.add_argument('-o', '--outdir', type=str, metavar='/path/to/output/directory/',
                                               help='Path to the directory where the information will be stored'
                                                    ' (default: ./)',
                                               required=False, default='.')

    parser_reads.set_defaults(func=reads_subcommand)
    parser_index.set_defaults(func=index_subcommand)
    parser_assembly.set_defaults(func=assembly_subcommand)
    parser_blast.set_defaults(func=blast_subcommand)

    return parser, parser_reads, parser_index, parser_assembly, parser_blast


def main():
    program_name = 'seq_typing.py'

    if sys.version_info[0] < 3:
        sys.exit('Must be using Python 3. Try calling "python3 {}"'.format(program_name))

    parser, _, _, _, _ = python_arguments(program_name, __version__)
    args = parser.parse_args()

    start_time = time.time()

    args.outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(args.outdir):
        os.makedirs(args.outdir)

    # Start logger
    logfile, time_str = utils.start_logger(args.outdir)

    script_path = utils.general_information(script_name=program_name, logfile=logfile, version=__version__,
                                            outdir=args.outdir, time_str=time_str)
    del script_path
    print('\n')

    folders_2_remove = []

    # Create modules pickles folder
    pickles_folder = os.path.join(args.outdir, 'pickles', '')
    if not os.path.isdir(pickles_folder):
        os.makedirs(pickles_folder)
    folders_2_remove.append(pickles_folder)

    # Run functions
    folders_2_remove_func, references_results, reference, references_headers = args.func(args)
    folders_2_remove.extend(folders_2_remove_func)

    # Parse results
    _, _, _, _, _ = parse_results.parse_results(references_results, reference, references_headers, args.outdir,
                                                args.minGeneCoverage, args.minDepthCoverage, args.typeSeparator)

    if not args.debug:
        for folder in folders_2_remove:
            utils.removeDirectory(folder)

    _ = utils.runTime(start_time)


if __name__ == "__main__":
    main()
