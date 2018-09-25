import os.path
import sys
from functools import partial

import modules.utils as utils


def check_db_exists(db_path):
    """
    Check if Blast DB exists

    Parameters
    ----------
    db_path : str
        Path to Blast DB files, e.g. /output/blast_db.nucl.sequences.fasta

    Returns
    -------
    db_exists : bool
        Tells if the Blast DB already exists
    original_file : bool
        Tells if the original fasta file from which the Blast DB was produced is present
    """
    db_exists = False
    files = [f for f in os.listdir(os.path.dirname(db_path))
             if not f.startswith('.')
             and os.path.isfile(os.path.join(os.path.dirname(db_path), f))]
    counter = 0
    original_file = False
    for file_found in files:
        if file_found == os.path.basename(db_path):
            original_file = True
        elif os.path.splitext(file_found)[0] == os.path.basename(db_path):
            counter += 1
    if counter > 0:
        db_exists = True
    return db_exists, original_file


def create_blast_db(db_sequences, db_output, db_type):
    """
    Creates a Blast DB

    Parameters
    ----------
    db_sequences : str
        Path to fasta file containing the sequences from which Blast DB will be created, e.g. /input/sequences.fasta
    db_output : str
        Path to Blast output DB files, e.g. /output/blast_db.nucl.sequences.fasta
    db_type : str
        Blast DB type. Can only be 'nucl' or 'prot'

    Returns
    -------
    run_successfully : bool
        Tells if the Blast DB was successfully created
    """

    run_successfully = False

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}). Use one of the following: nucl, prot'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(db_output)):
        os.makedirs(os.path.dirname(db_output))
    else:
        db_exists, original_file = check_db_exists(db_output)
        if db_exists and original_file:
            print('Blast DB already found at {db_output}'
                  ' for {db_sequences}'.format(db_output=os.path.dirname(db_output),
                                               db_sequences=db_sequences,
                                               file_found=os.path.basename(db_sequences)))

            run_successfully = True
        elif db_exists and not original_file:
            print('The original fasta file ({db_sequences}) from which the Blast DB was produced is not present. Make'
                  ' sure it is found in {db_dir} and it is'
                  ' named {original_file_name}'.format(db_dir=os.path.dirname(db_output), db_sequences=db_sequences,
                                                       original_file_name=os.path.basename(db_output)))
        elif not db_exists and original_file:
            print('The original fasta file ({db_sequences}) from which the Blast DB was supposed to be produced is'
                  ' present in Blast DB directory ({db_dir}), but no Blast DB files were found'
                  ' there.'.format(db_sequences=db_sequences, db_dir=os.path.dirname(db_output)))
        else:
            run_successfully, _, _ = utils.runCommandPopenCommunicate(['makeblastdb', '-parse_seqids', '-dbtype',
                                                                       db_type, '-in', db_sequences, '-out', db_output],
                                                                      False, None, True)
            if run_successfully:
                from shutil import copyfile
                copyfile(db_sequences, db_output)

    return run_successfully


def run_blast_command(query_file, blast_db, db_type, blast_output, threads=1):
    """
    Run Blast: blastn or blastp

    Parameters
    ----------
    query_file : str
        Path to fasta file containing the query sequences, e.g. /input/queries.fasta
    blast_db : str
        Path to Blast DB files, e.g. /input/blast_db.nucl.sequences.fasta
    db_type : str
        Blast DB type. Can only be 'nucl' or 'prot'
    blast_output : str
        Path to Blast output tabular file
    threads : int
        Number of CPUs to use during Blast

    Returns
    -------
    run_successfully : bool
        Tells if the Blast DB was successfully created
    """

    # Check Blast DB type
    if db_type not in ('nucl', 'prot'):
        exit('Wrong Blast DB type provided ({db_type}).'
             ' Use one of the following: "nucl", "prot"'.format(db_type=db_type))

    # Check if db_output dir exists
    if not os.path.isdir(os.path.dirname(blast_output)):
        os.makedirs(os.path.dirname(blast_output))

    command = ['', '-query', query_file, '-db', blast_db, '-out', blast_output, '-outfmt', '', '-dust', 'no',
               '-culling_limit', '1', '-num_threads', str(threads)]

    if db_type == 'nucl':
        command[0] = 'blastn -task blastn'
    else:
        command[0] = 'blastp -task blastp'

    command[8] = "'7 qseqid qlen sseqid slen qstart qend sstart send evalue length pident nident mismatch gaps'"

    run_successfully, _, _ = utils.runCommandPopenCommunicate(command, False, None, True)

    return run_successfully


def parse_blast_output(blast_output):
    """
    Parse Blast output

    Parameters
    ----------
    blast_output : str
        Path to Blast output tabular file

    Returns
    -------
    output_blast : dict
        Dictionary with Blast results cleaned and almost already formated for parse_results.py: subject as key and
        values similar to ReMatCh results
    """

    # Blast fields
    # 'query_id', 'query_length', 'subject_id', 'subject_length', 'q_start', 'q_end', 's_start', 's_end', 'evalue',
    # 'alignment_length', 'percentage_identity', 'identical', 'mismatches', 'gaps'

    output_blast = {}

    # Read Blast output
    with open(blast_output, 'rt') as reader:
        for line in reader:
            line = line.rstrip('\r\n')
            if len(line) > 0:
                if not line.startswith('#'):
                    line = line.split('\t')

                    # By subject
                    if line[2] not in output_blast:
                        output_blast[line[2]] = {'header': line[2],
                                                 'gene_coverage': int(line[9]) / float(line[3]) * 100,
                                                 'gene_low_coverage': 0,
                                                 'gene_number_positions_multiple_alleles': 0,
                                                 'gene_mean_read_coverage': 1,
                                                 'gene_identity': float(line[10]),
                                                 'q_start': int(line[4]),
                                                 'q_end': int(line[5]), 's_start': int(line[6]),
                                                 's_end': int(line[7]), 'evalue': float(line[8]),
                                                 'gaps': int(line[13]), 'query': line[0],
                                                 'alignment_length': int(line[9]), 'subject_length': int(line[3])}
                    else:
                        previous_blast_results = output_blast[line[2]]
                        present_blast_results = {'header': line[2],
                                                 'gene_coverage': int(line[9]) / float(line[3] * 100),
                                                 'gene_low_coverage': 0, 'gene_number_positions_multiple_alleles': 0,
                                                 'gene_mean_read_coverage': 1, 'gene_identity': float(line[10]),
                                                 'q_start': int(line[4]), 'q_end': int(line[5]),
                                                 's_start': int(line[6]), 's_end': int(line[7]),
                                                 'evalue': float(line[8]),
                                                 'gaps': int(line[13]), 'query': line[0],
                                                 'alignment_length': int(line[9]), 'subject_length': int(line[3])}

                        to_change = False
                        if previous_blast_results['alignment_length'] < int(line[9]):
                            to_change = True
                        # elif previous_blast_results['alignment_length'] == int(line[9]) and \
                        #         previous_blast_results['gene_coverage'] <= present_blast_results['gene_coverage']:
                        #     if previous_blast_results['gene_coverage'] < present_blast_results['gene_coverage']:
                        #         to_change = True
                        elif previous_blast_results['alignment_length'] == int(line[9]) and \
                                previous_blast_results['gene_identity'] <= float(line[10]):
                            if previous_blast_results['gene_identity'] < float(line[10]):
                                to_change = True
                            elif previous_blast_results['gene_identity'] == float(line[10]) and \
                                    previous_blast_results['evalue'] >= float(line[8]):
                                if previous_blast_results['evalue'] > float(line[8]):
                                    to_change = True
                                elif previous_blast_results['evalue'] == float(line[8]) and \
                                        previous_blast_results['gaps'] > int(line[13]):
                                    to_change = True

                        if to_change:
                            output_blast[line[2]] = present_blast_results

    return output_blast


module_timer = partial(utils.timer, name='Module Blast')


@module_timer
def run_blast(blast_db_path, outdir, blast_type, query_fasta_file):
    """
    Parse Blast output

    Parameters
    ----------
    blast_db_path : str
        Path to Blast DB files, e.g. /input/blast_db.nucl.sequences.fasta. If Blast DB is not found, it will be created
        under the outdir folder
    outdir : str
        Path to output directory
    blast_type : str
        Blast DB type. Can only be 'nucl' or 'prot'
    query_fasta_file : str
        Path to fasta file containing the query sequences, e.g. /input/queries.fasta

    Returns
    -------
    folders_2_remove : list
        List of folders that can be removed at the end
    blast_results : dict
        Dictionary with Blast results cleaned and almost already formated for parse_results.py: subject as key and
        values similar to ReMatCh results
    blast_db_path : str
        Path to Blast DB used
    headers_correspondence : dict
        Dictionary sequence headers correspondence if headers were changed. Keys are the new headers and values the
        original ones. If headers were not changed, keys and vules are the same.
    """

    folders_2_remove = []

    # Check Blast DB
    db_exists, original_file = check_db_exists(blast_db_path)
    if not db_exists:
        blast_db = os.path.join(outdir, 'blast_db', '{blast_DB}'.format(blast_DB=os.path.basename(blast_db_path)))
        folders_2_remove.append(os.path.dirname(blast_db))

        if not os.path.isdir(os.path.dirname(blast_db)):
            os.makedirs(os.path.dirname(blast_db))

        db_exists = create_blast_db(blast_db_path, blast_db, blast_type)
        if db_exists:
            blast_db_path = str(blast_db)
            original_file = True
    elif db_exists and not original_file:
        sys.exit('Original fasta file from which the Blast DB was produced ({blast_db_path}) is missing'
                 ' from {blast_db_dir}'.format(blast_db_path=os.path.basename(blast_db_path),
                                               blast_db_dir=os.path.dirname(blast_db_path)))

    # Run Blast
    if db_exists and original_file:
        sys.path.append('..')
        from seq_typing import parse_reference
        _, headers_correspondence = parse_reference(blast_db_path, [])

        blast_output = os.path.join(outdir, 'blast_out', 'results.tab')
        folders_2_remove.append(os.path.dirname(blast_output))

        if not os.path.isdir(os.path.dirname(blast_output)):
            os.makedirs(os.path.dirname(blast_output))
        run_successfully = run_blast_command(query_fasta_file, blast_db_path, blast_type, blast_output)

        if run_successfully:
            blast_results = parse_blast_output(blast_output)
        else:
            sys.exit('Blast was not run successfully')
    else:
        sys.exit('It was not found any Blast DB and/or the original fasta file from which the Blast DB was produced')

    return folders_2_remove, blast_results, blast_db_path, headers_correspondence
