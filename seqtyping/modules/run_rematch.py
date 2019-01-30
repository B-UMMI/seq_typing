from functools import partial
import os
import sys
from itertools import groupby as itertools_groupby

try:
    import modules.utils as utils
except ImportError:
    from seqtyping.modules import utils as utils


def remove_alignment(alignment_dir):
    if os.path.isdir(alignment_dir):
        files = [f for f in os.listdir(alignment_dir) if
                 not f.startswith('.') and
                 os.path.isfile(os.path.join(alignment_dir, f)) and
                 f.endswith(('.bam', '.sam', '.cram'))]
        for file_found in files:
            file_found = os.path.join(alignment_dir, file_found)
            os.remove(file_found)


def remove_reference_stuff(outdir, reference_file):
    files = [f for f in os.listdir(outdir) if not f.startswith('.') and os.path.isfile(os.path.join(outdir, f))]
    for file_found in files:
        if file_found.startswith(os.path.splitext(os.path.basename(reference_file))[0]):
            file_found = os.path.join(outdir, file_found)
            os.remove(file_found)


def clean_rematch_folder(consensus_files, reference_file, outdir, doNotRemoveConsensus, debug_mode_true):
    if not debug_mode_true:
        if consensus_files is not None:
            if not doNotRemoveConsensus:
                for consensus_type, file_path in consensus_files.items():
                    if os.path.isfile(file_path):
                        os.remove(file_path)
        remove_alignment(os.path.join(outdir, 'rematch_module', ''))
        remove_reference_stuff(outdir, reference_file)


problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]


def clean_header(header):
    new_header = str(header)
    if any(x in header for x in problematic_characters):
            for x in problematic_characters:
                new_header = new_header.replace(x, '_')
    return header, new_header


def get_sequence_information(fasta_file):
    headers = {}
    sequence_dict = {}
    headers_changed = False

    sequence_counter = 0

    reader = open(fasta_file, mode='rt', newline=None)
    fasta_iter = (g for k, g in itertools_groupby(reader, lambda x: x.startswith('>')))
    for header in fasta_iter:
        original_header, new_header = clean_header(header.__next__()[1:].rstrip('\r\n'))
        if new_header in headers:
            sys.exit('Found duplicated sequence'
                     ' headers: {original_header}'.format(original_header=original_header))
        seq = ''.join(s.rstrip('\r\n') for s in fasta_iter.__next__())
        sequence_counter += 1
        sequence_dict[sequence_counter] = {'header': new_header, 'sequence': seq, 'length': len(seq)}
        headers[new_header] = str(original_header)
        if new_header != original_header:
            headers_changed = True
    reader.close()

    return sequence_dict, headers, headers_changed


def clean_headers_reference_file(reference_file, outdir):
    new_reference_file = reference_file
    sequences, headers, headers_changed = get_sequence_information(reference_file)
    if headers_changed:
        utils.Bcolors_print(str('At least one of the following characters was found in sequences headers: {}.\n'
                                'Replacing those with _'.format(problematic_characters) + '\n'), 'UNDERLINE')
        new_reference_file = \
            os.path.join(outdir, os.path.splitext(os.path.basename(reference_file))[0] + '.headers_renamed.fasta')
        with open(new_reference_file, 'wt') as writer:
            for i in sequences:
                writer.write('>' + sequences[i]['header'] + '\n')
                fasta_sequence_lines = utils.chunkstring(sequences[i]['sequence'], 80)
                for line in fasta_sequence_lines:
                    writer.write(line + '\n')
    return new_reference_file, headers, sequences


def rematch_for_different_references(fastq, references_files, threads, outdir, extraSeq, minCovPresence, minCovCall,
                                     minFrequencyDominantAllele, minGeneCoverage, debug, minGeneIdentity,
                                     rematch_module, doNotRemoveConsensus, bowtie_algorithm='--very-sensitive-local',
                                     max_number_mapped_location=1, clean_run_rematch=False, save_new_allele=False):
    if save_new_allele:
        os.makedirs(os.path.join(outdir, 'new_allele', ''))

    references_results = {}
    for x, reference in enumerate(references_files):
        reference_name = os.path.basename(reference) + '_' + str(x)
        ref_dir = os.path.join(outdir, reference_name, '')
        os.makedirs(ref_dir)
        header_gene_list, seq_reference_dict = utils.extractVariableFromPickle(reference + '.pkl')

        not_write_consensus_rematch = not doNotRemoveConsensus and not save_new_allele

        time_taken, run_successfully, data_by_gene, sample_data_general, consensus_files, consensus_sequences = \
            rematch_module.run_rematch_module('sample', fastq, reference, threads, ref_dir, extraSeq,
                                              minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage,
                                              debug, max_number_mapped_location, minGeneIdentity, 'first', 7, 'none',
                                              seq_reference_dict, 'X', bowtie_algorithm, None, header_gene_list,
                                              not_write_consensus_rematch, clean_run=clean_run_rematch)
        if run_successfully:
            pickleFile = os.path.join(outdir, str(reference_name + '.pkl'))
            utils.saveVariableToPickle(data_by_gene, pickleFile)
            references_results[reference] = pickleFile
        else:
            sys.exit('Something went wrong while running ReMatCh for reference {reference}'.format(reference=reference))

        if save_new_allele:
            from shutil import copyfile
            copyfile(os.path.join(ref_dir, 'sample.noMatter.fasta'), os.path.join(outdir, 'new_allele', reference))

        clean_rematch_folder(consensus_files, reference, ref_dir, doNotRemoveConsensus, debug)
        if not debug and not doNotRemoveConsensus:
            utils.removeDirectory(ref_dir)
    return references_results


def run_bowtie_inspect(index_without_sufix, outdir):
    """
    Create the reference fasta file from Bowtie2 index files

    Parameters
    ----------
    index_without_sufix : str
        Path to the basename of the index. The basename is name of any of the index files but with the .X.bt2 or
        .rev.X.bt2 suffix omitted.
    outdir : str
        Path to output directory where the fasta file will be stored

    Returns
    -------
    run_successfully : bool
        Tells if the reference fasta file was successfully created
    reference_file : str
        Path to the newly created reference fasta file
    """

    reference_file = os.path.join(outdir, os.path.basename(index_without_sufix))
    run_successfully, _, _ = utils.runCommandPopenCommunicate(['bowtie2-inspect', '--across', 80, index_without_sufix,
                                                               '>', reference_file],
                                                              True, None, True)
    return run_successfully, reference_file


def run_bowtie_build(reference_file, outdir, threads=1):
    """
    Create Bowtie2 index files

    Parameters
    ----------
    reference_file : str
        Path to the reference sequences file
    outdir : str
        Path to output directory where the fasta file will be stored
    threads : int, default 1
        Number of CPUs to use while creating Bowtie2 index
    Returns
    -------
    run_successfully : bool
        Tells if the reference fasta file was successfully created
    """

    index_basename = os.path.join(outdir, os.path.basename(reference_file))
    run_successfully, _, _ = utils.runCommandPopenCommunicate(['bowtie2-build', '--threads', str(threads),
                                                               reference_file, index_basename],
                                                              False, None, True)
    return run_successfully


module_timer = partial(utils.timer, name='Module ReMatCh')


@module_timer
def run_rematch(rematch_script, outdir, references_files, fastq, threads, extraSeq, minCovPresence, minCovCall,
                minFrequencyDominantAllele, minGeneCoverage, minGeneIdentity, debug, doNotRemoveConsensus,
                bowtie_algorithm='--very-sensitive-local', max_number_mapped_location=1, clean_run_rematch=False,
                save_new_allele=False):
    module_dir = os.path.join(outdir, 'rematch', '')
    utils.removeDirectory(module_dir)
    os.makedirs(module_dir)

    sys.path.append(os.path.join(os.path.dirname(rematch_script), 'modules'))
    import rematch_module

    references_results = rematch_for_different_references(fastq, references_files, threads, module_dir, extraSeq,
                                                          minCovPresence, minCovCall, minFrequencyDominantAllele,
                                                          minGeneCoverage, debug, minGeneIdentity, rematch_module,
                                                          doNotRemoveConsensus, bowtie_algorithm=bowtie_algorithm,
                                                          max_number_mapped_location=max_number_mapped_location,
                                                          clean_run_rematch=clean_run_rematch,
                                                          save_new_allele=save_new_allele)

    return references_results, module_dir
