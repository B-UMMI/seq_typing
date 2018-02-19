import functools
import os
import sys

import modules.utils as utils


def remove_alignment(alignment_dir):
    if os.path.isdir(alignment_dir):
        files = [f for f in os.listdir(alignment_dir) if not f.startswith('.') and os.path.isfile(os.path.join(alignment_dir, f)) and f.endswith(('.bam', '.sam', '.cram'))]
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


def clean_headers_reference_file(reference_file, outdir, extraSeq, rematch_module):
    problematic_characters = ["|", " ", ",", ".", "(", ")", "'", "/", ":"]
    print('Checking if reference sequences contain ' + str(problematic_characters) + '\n')
    headers_changed = False
    new_reference_file = reference_file
    sequences, headers, headers_changed = rematch_module.get_sequence_information(reference_file, extraSeq)
    for i in sequences:
        if any(x in sequences[i]['header'] for x in problematic_characters):
            for x in problematic_characters:
                sequences[i]['header'] = sequences[i]['header'].replace(x, '_')
            headers_changed = True
    if headers_changed:
        utils.Bcolors_print(str('At least one of the those characters was found. Replacing those with _' + '\n'), 'UNDERLINE')
        new_reference_file = os.path.join(outdir, os.path.splitext(os.path.basename(reference_file))[0] + '.headers_renamed.fasta')
        with open(new_reference_file, 'wt') as writer:
            for i in sequences:
                writer.write('>' + sequences[i]['header'] + '\n')
                fasta_sequence_lines = rematch_module.chunkstring(sequences[i]['sequence'], 80)
                for line in fasta_sequence_lines:
                    writer.write(line + '\n')
    return new_reference_file, headers, sequences


def rematch_for_different_references(fastq, references_files, threads, outdir, extraSeq, minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage, debug, minGeneIdentity, rematch_module, doNotRemoveConsensus):
    references_results = {}
    for x, reference in enumerate(references_files):
        reference_name = os.path.basename(reference) + '_' + str(x)
        ref_dir = os.path.join(outdir, reference_name, '')
        os.makedirs(ref_dir)
        reference_file, gene_list_reference, reference_dict = clean_headers_reference_file(reference, ref_dir, extraSeq, rematch_module)
        time_taken, run_successfully, data_by_gene, sample_data_general, consensus_files, consensus_sequences = rematch_module.runRematchModule('sample', fastq, reference_file, threads, ref_dir, extraSeq, minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage, False, debug, 1, minGeneIdentity, 'first', 7, 'none', reference_dict, 'X', None, gene_list_reference, not doNotRemoveConsensus)
        if run_successfully:
            pickleFile = os.path.join(outdir, str(reference_name + '.pkl'))
            utils.saveVariableToPickle(data_by_gene, pickleFile)
            references_results[reference] = pickleFile
        else:
            sys.exit('Something went wrong while running ReMatCh for reference {reference}'.format(reference=reference))
        clean_rematch_folder(consensus_files, reference_file, ref_dir, doNotRemoveConsensus, debug)
        if not debug and not doNotRemoveConsensus:
            utils.removeDirectory(ref_dir)
    return references_results


module_timer = functools.partial(utils.timer, name='Module ReMatCh')


@module_timer
def run_rematch(rematch_script, outdir, references_files, fastq, threads, extraSeq, minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage, minGeneIdentity, debug, doNotRemoveConsensus):
    module_dir = os.path.join(outdir, 'rematch', '')
    utils.removeDirectory(module_dir)
    os.makedirs(module_dir)

    sys.path.append(os.path.join(os.path.dirname(rematch_script), 'modules', ''))
    from past import autotranslate
    autotranslate(['rematch_module', 'utils'])
    import rematch_module

    references_results = rematch_for_different_references(fastq, references_files, threads, module_dir, extraSeq, minCovPresence, minCovCall, minFrequencyDominantAllele, minGeneCoverage, debug, minGeneIdentity, rematch_module, doNotRemoveConsensus)

    return references_results, module_dir
