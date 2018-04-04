# seq_typing

Determine which reference sequence is more likely to be present in a given sample

---

* [Rational](#rational)
* [Requirements](#requirements)
* [Dependencies](#dependencies)
  * [Install dependencies](#install-dependencies)
* [Install seq_typing](#install-seq_typing)
* [Usage](#usage)
  * [Species serotyping](#species-serotyping)
  * [Usage examples](#usage-examples)
* [Outputs](#outputs)
* [Contact](#contact)

## Rational

**seq_typing** is a software to determine a given sample type using a read mapping approach against a set of reference sequences. Sample's reads are mapped to the given reference sequences and, based on the length of the sequence covered and it's depth of coverage, **seq_typing** decides which reference sequence is more likely to be present and returns the type associated with such sequence.

## Requirements

* Illumina Fastq files

## Dependencies

* Python 3
* Python module _**future**_
* [ReMatCh](https://github.com/B-UMMI/ReMatCh)

### Install dependencies

Python using [Conda](https://conda.io/) (Python 3):

````bash
conda create --name seq_typing python=3 future
````

ReMatCh:
````bash
git clone https://github.com/B-UMMI/ReMatCh.git
cd ReMatCh

# Temporarily add ReMatCh to the PATH
export PATH="$(pwd -P):$PATH"

# Permanently add ReMatCh to the PATH
echo export PATH="$(pwd -P):$PATH" >> ~/.profile
````

## Install seq_typing

````bash
git clone https://github.com/B-UMMI/seq_typing.git
cd seq_typing

# Temporarily add seq_typing to the PATH
export PATH="$(pwd -P):$PATH"

# Permanently add seq_typing to the PATH
echo export PATH="$(pwd -P):$PATH" >> ~/.profile
````

## Usage

````
usage: seq_typing.py [-h] [--version] -f /path/to/input/file.fq.gz
                     [/path/to/input/file.fq.gz ...]
                     [-r /path/to/reference_sequence.fasta [/path/to/reference_sequence.fasta ...]
                     | -s Escherichia coli] [-o /path/to/output/directory/]
                     [-j N] [--mapRefTogether] [--typeSeparator _]
                     [--extraSeq N] [--minCovPresence N] [--minCovCall N]
                     [--minGeneCoverage N] [--doNotRemoveConsensus] [--debug]
                     [--beginning] [--notClean]

Determine which reference sequence is more likely to be present in a given
sample

optional arguments:
  -h, --help            show this help message and exit
  --version             Version information
  -r /path/to/reference_sequence.fasta [/path/to/reference_sequence.fasta ...], --reference /path/to/reference_sequence.fasta
                        Fasta file containing reference sequences. If more
                        than one file is passed, a reference sequence for each
                        file will be determined. Give the files name in the
                        same order that the type must be determined. (default:
                        None)
  -s Escherichia coli, --species Escherichia coli
                        Name of the species with reference sequences provided
                        together with seq_typing.py for serotyping (default:
                        None)

Required options:
  -f /path/to/input/file.fq.gz [/path/to/input/file.fq.gz ...], --fastq /path/to/input/file.fq.gz
                        Path to single OR paired-end fastq files. If two files
                        are passed, they will be assumed as being the paired
                        fastq files (default: None)

General facultative options:
  -o /path/to/output/directory/, --outdir /path/to/output/directory/
                        Path to the directory where the information will be
                        stored (default: .)
  -j N, --threads N     Number of threads to use (default: 1)
  --mapRefTogether      Map the reads against all references together
                        (default: False)
  --typeSeparator _     Last single character separating the general sequence
                        header from the last part containing the type
                        (default: _)
  --extraSeq N          Sequence length added to both ends of target sequences
                        (usefull to improve reads mapping to the target one)
                        that will be trimmed in ReMatCh outputs (default: 0)
  --minCovPresence N    Reference position minimum coverage depth to consider
                        the position to be present in the sample (default: 5)
  --minCovCall N        Reference position minimum coverage depth to perform a
                        base call (default: 10)
  --minGeneCoverage N   Minimum percentage of target reference sequence
                        covered to consider a sequence to be present (value
                        between [0, 100]) (default: 60)
  --doNotRemoveConsensus
                        Do not remove ReMatCh consensus sequences (default:
                        False)
  --debug               DeBug Mode: do not remove temporary files (default:
                        False)
  --beginning           Start seq_typing.py from the beggining (default:
                        False)
  --notClean            Do not remove intermediate files (default: False)
````

### Species serotyping

For the following species, references sequences are provided for serotyping.
* Escherichia coli
* Haemophilus influenzae
* Streptococcus agalactiae

Use `--species` option with one of those species

### Usage examples

Serotyping _Haemophilus influenzae_ using provided references sequences (that uses only one reference sequences file):
````bash
# Activate Conda environment (when using Python via Conda)
source activate seq_typing

seq_typing.py --species Haemophilus influenzae \
              --fastq sample_1.fq.gz sample_2.fq.gz \
              --outdir sample_out/ \
              --threads 2
````

Serotyping _Escherichia coli_ using provided references sequences (that uses two reference sequences files):
````bash
# Activate Conda environment (when using Python via Conda)
source activate seq_typing

seq_typing.py --species Escherichia coli \
              --fastq sample_1.fq.gz sample_2.fq.gz \
              --outdir sample_out/ \
              --threads 2 \
              --mapRefTogether
````

Type one sample with a given set of references sequences:
````bash
# Activate Conda environment (when using Python via Conda)
source activate seq_typing

seq_typing.py --reference references/Ecoli/O_type.fasta references/Ecoli/H_type.fasta \
              --fastq sample_1.fq.gz sample_2.fq.gz \
              --outdir sample_out/ \
              --threads 2 \
              --mapRefTogether
````

## Outputs

__seq_typing.report.txt__  
Text file with the typing result. If it was not possible to determine a type for a given reference file, `NT` (for None Typeable) will be returned for that file.

Example of _E. coli_ serotyping:  
`O157:H7`

__seq_typing.report_types.tab__  
Tabular file with detailed results:
* _sequence_type_: type of the results reported. Three values can be found here. `selected` for the reference sequences selected for the reported typing result. `other_probable_type` for other reference sequences that could have been selected because fulfill selection thresholds. `most_likely` for the most likely reference sequences when no reference sequences fulfill selection thresholds.
* _reference_file_: the reference file where the sequences came from.
* _sequence_: reference sequences name
* _sequenced_covered_: percentage of reference sequences covered
* _coverage_depth_: mean reference sequences depth of coverage of the positions present
* _sequence_identity_: percentage identity of reference sequences covered

Example of _E. coli_ serotyping:

| #sequence_type      | reference_file | sequence              | sequenced_covered | coverage_depth     | sequence_identity |
|---------------------|----------------|-----------------------|-------------------|--------------------|-------------------|
| selected            | O_type.fasta   | wzy_192_AF529080_O26  | 100.0             | 281.95405669599216 | 100.0             |
| selected            | H_type.fasta   | fliC_269_AY337465_H11 | 99.4546693933197  | 51.76490747087046  | 99.86291980808772 |
| other_probable_type | O_type.fasta   | wzx_208_AF529080_O26  | 100.0             | 223.3072050673001  | 100.0             |
| other_probable_type | H_type.fasta   | fliC_276_AY337472_H11 | 98.84117246080436 | 37.52551724137931  | 99.86206896551724 |

__run.*.log__  
Running log file.  

## Contact

Miguel Machado  
<mpmachado@medicina.ulisboa.pt>
