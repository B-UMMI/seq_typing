# seq_typing

Determines which reference sequence is more likely to be present in a given sample

---

* [Rational](#rational)
* [Requirements](#requirements)
* [Dependencies](#dependencies)
  * [Install dependencies](#install-dependencies)
* [Install seq_typing](#install-seq_typing)
* [Usage](#usage)
  * [General use](#general-use)
    * [General info](#general-info)
      * [_reads_ module](#reads-module)
      * [_blast_ module](#blast-module)
      * [_assembly_ module](#assembly-module)
    * [Organisms typing](#organisms-typing)
    * [Usage examples](#usage-examples)
      * [Reads](#reads)
      * [Assemblies](#assemblies)
  * [E. coli stx subtyping](#e-coli-stx-subtyping)
    * [General usage](#general-usage)
    * [Update stx references](#update-stx-references)
* [Outputs](#outputs)
* [Contact](#contact)

## Rational

**seq_typing** is a software to determine a given sample type using a read mapping approach or sequence Blast search against a set of reference sequences. The sample's reads are mapped to the given reference sequences and, based on the length of the sequence covered and it's depth of coverage, **seq_typing** decides which reference sequence is the most likely to be present, and returns the type associated with such sequence. When using sequences fasta files, similar decision rules are applied, but results are first cleaned to get the best hit for each DB sequence (based on alignment length, similarity, E-value and number of gaps).

## Requirements

* Illumina Fastq files  
_OR_
* Sequence fasta file

## Dependencies

* Python 3
* Fastq files:
  * Python module _**future**_
  * [ReMatCh](https://github.com/B-UMMI/ReMatCh)
* Fasta file:
  * [Blast+](https://www.ncbi.nlm.nih.gov/books/NBK279690/)

### Install dependencies

Python using [Conda](https://conda.io/) (Python 3, with _future_ module):

```bash
conda create --name seq_typing python=3 future
```

ReMatCh:
```bash
git clone https://github.com/B-UMMI/ReMatCh.git
cd ReMatCh

# Temporarily add ReMatCh to the PATH
export PATH="$(pwd -P):$PATH"

# Permanently add ReMatCh to the PATH
echo export PATH="$(pwd -P):$PATH" >> ~/.profile
```

Blast+:
```bash
wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/ncbi-blast-*-x64-linux.tar.gz
tar xf ncbi-blast-*-x64-linux.tar.gz
rm ncbi-blast-*-x64-linux.tar.gz
cd ncbi-blast-*/bin

# Temporarily add Blast binaries to the PATH
export PATH="$(pwd -P):$PATH"

# Permanently add Blast binaries to the PATH
echo export PATH="$(pwd -P):$PATH" >> ~/.profile
```

## Install seq_typing

```bash
git clone https://github.com/B-UMMI/seq_typing.git
cd seq_typing

# Temporarily add seq_typing to the PATH
export PATH="$(pwd -P):$PATH"

# Permanently add seq_typing to the PATH
echo export PATH="$(pwd -P):$PATH" >> ~/.profile
```

## Usage

### General use

#### General info
```
usage: seq_typing.py [-h] [--version] {reads,assembly,blast} ...

Determines which reference sequence is more likely to be present in a given
sample

optional arguments:
  -h, --help            Show this help message and exit
  --version             Version information

Subcommands:
  Valid subcommands

  {reads,assembly,blast}
                        Additional help
    reads               reads --help
    assembly            assembly --help
    blast               blast --help
```
* [_reads_ module](#reads-module):  
  Run seq_typing.py using fastq files.
* [_blast_ module](#blast-module):  
  Creates Blast DB. This is useful when running the same DB sequence file for different samples.
* [_assembly_ module](#assembly-module):  
  Run seq_typing.py using a fasta file. If running multiple samples using the same DB sequence file, consider use first _seq_typing.py blast_ module.

##### _reads_ module  
Run seq_typing.py using fastq files.
```
usage: seq_typing.py reads [-h]
                           -f /path/to/input/file.fq.gz ...
                           -r /path/to/reference_sequence.fasta ... | --org escherichia coli
                           [-o /path/to/output/directory/] [-j N]
                           [--mapRefTogether] [--typeSeparator _]
                           [--extraSeq N] [--minCovPresence N]
                           [--minCovCall N] [--minGeneCoverage N]
                           [--minDepthCoverage N] [--minGeneIdentity N]
                           [--doNotRemoveConsensus] [--debug] [--resume]
                           [--notClean]

Run seq_typing.py using fastq files

optional arguments:
  -h, --help            show this help message and exit

Required options:
  -f --fastq /path/to/input/file.fq.gz ...
                        Path to single OR paired-end fastq files. If two files
                        are passed, they will be assumed as being the paired
                        fastq files

Required one of the following options:
  -r --reference /path/to/reference_sequence.fasta ...
                        Fasta file containing reference sequences. If more
                        than one file is passed, a type for each file will be
                        determined. Give the files name in the same order that
                        the type must be determined.
  --org escherichia coli
                        Name of the organism with reference sequences provided
                        together with seq_typing.py for typing
                        ("reference_sequences" folder)

General facultative options:
  -o --outdir /path/to/output/directory/
                        Path to the directory where the information will be
                        stored (default: ./)
  -j --threads N        Number of threads to use (default: 1)
  --mapRefTogether      Map the reads against all references together
  --typeSeparator _     Last single character separating the general sequence
                        header from the last part containing the type (default: _)
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
  --minDepthCoverage N  Minimum depth of coverage of target reference sequence
                        to consider a sequence to be present (default: 2)
  --minGeneIdentity N   Minimum percentage of identity of reference sequence
                        covered to consider a gene to be present (value
                        between [0, 100]). One INDEL will be considered as one
                        difference (default: 80)
  --doNotRemoveConsensus
                        Do not remove ReMatCh consensus sequences
  --debug               Debug mode: do not remove temporary files
  --resume              Resume seq_typing.py reads
```

##### _blast_ module  
Creates Blast DB.  
This is useful when running the same DB sequence file for different samples.
```
usage: seq_typing.py blast [-h]
                           -t nucl
                           -f /path/to/db.sequences.fasta ... | --org escherichia coli
                           [-o /path/to/output/directory/]

Creates Blast DB. This is useful when running the same DB sequence file for
different samples.

optional arguments:
  -h, --help            show this help message and exit

Required one of the following options:
  -f --fasta /path/to/db.sequences.fasta ...
                        Path to DB sequence file. If more than one file is
                        passed, a Blast DB for each file will be created.
  --org escherichia coli
                        Name of the organism with DB sequence file provided
                        ("reference_sequences" folder) together with
                        seq_typing.py for typing

Required option for --fasta:
  -t nucl, --type nucl  Blast DB type (available options: nucl, prot)

General facultative options:
  -o --outdir /path/to/output/directory/
                        Path to the directory where the information will be
                        stored (default: ./)
```

##### _assembly_ module  
Run **seq_typing** using a fasta file.  
If running multiple samples using the same DB sequence file, consider use first _seq_typing.py blast_ module.
```
usage: seq_typing.py assembly [-h]
                              -f /path/to/query/assembly_file.fasta
                              -b /path/to/Blast/db.sequences.file ... -t nucl | --org escherichia coli
                              [-o /path/to/output/directory/] [-j N]
                              [--typeSeparator _] [--minGeneCoverage N]
                              [--minGeneIdentity N] [--debug]

Run seq_typing.py using a fasta file. If running multiple samples using the
same DB sequence file, consider use first "seq_typing.py blast"
module.

optional arguments:
  -h, --help            show this help message and exit

Required options:
  -f /path/to/query/assembly_file.fasta, --fasta /path/to/query/assembly_file.fasta
                        Path to fasta file containing the query sequences from
                        which the types should be assessed

Required one of the following options:
  -b --blast /path/to/Blast/db.sequences.file ...
                        Path to DB sequence file. If Blast DB was already
                        produced only provide the file that do not end with
                        ".n*" something (do not use for example
                        /blast_db.sequences.fasta.nhr). If no Blast DB is
                        found for the DB sequence file, one will be created in
                        --outdir. If more than one Blast DB file is passed, a
                        type for each file will be determined. Give the files
                        in the same order that the type must be determined.
  --org escherichia coli
                        Name of the organism with DB sequence file provided
                        ("reference_sequences" folder) together with
                        seq_typing.py for typing

Required option for --blast:
  -t --type nucl        Blast DB type (available options: nucl, prot)

General facultative options:
  -o --outdir /path/to/output/directory/
                        Path to the directory where the information will be
                        stored (default: ./)
  -j --threads N        Number of threads to use (default: 1)
  --typeSeparator _     Last single character separating the general sequence
                        header from the last part containing the type (default: _)
  --minGeneCoverage N   Minimum percentage of target reference sequence
                        covered to consider a sequence to be present (value
                        between [0, 100]) (default: 60)
  --minGeneIdentity N   Minimum percentage of identity of reference sequence
                        covered to consider a gene to be present (value
                        between [0, 100]) (default: 80)
  --debug               Debug mode: do not remove temporary files
```

#### Organisms typing

For the following organisms, references sequences are provided for serotyping.
* _Escherichia coli_
* _Haemophilus influenzae_
* _Streptococcus agalactiae_
* Dengue virus (with genotype information)

Use `--org` option with one of those organisms

#### Usage examples

##### Reads

Serotyping _Haemophilus influenzae_ using provided references sequences (that uses only one reference sequences file):
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

seq_typing.py reads --org Haemophilus influenzae \
                    --fastq sample_1.fq.gz sample_2.fq.gz \
                    --outdir sample_out/ \
                    --threads 2
```

Serotyping _Escherichia coli_ using provided references sequences (that uses two reference sequences files):
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

seq_typing.py reads --org Escherichia coli \
                    --fastq sample_1.fq.gz sample_2.fq.gz \
                    --outdir sample_out/ \
                    --threads 2 \
                    --mapRefTogether
```

Type one sample with a users own set of references sequences (using for example single-end reads):
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

seq_typing.py reads --reference references/Ecoli/O_type.fasta references/Ecoli/H_type.fasta \
                    --fastq sample.fq.gz \
                    --outdir sample_out/ \
                    --threads 2 \
                    --mapRefTogether
```

##### Assemblies

Type _Dengue virus_ using assemblies with provided reference sequences (uses only one reference sequences file):
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

seq_typing.py assembly --org Dengue virus \
                       --fasta sample.fasta \
                       --outdir sample_out/ \
                       --threads 2
```

When running the same database for different samples, a single Blast database should be produce first to speed up the analysis.  
Example using _Escherichia coli_ provided reference sequences (that uses two reference sequences files):
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

seq_typing.py blast --org Escherichia coli \
                    --outdir db_out/

# Run seq_typing using created database
seq_typing.py assembly --blast db_out/1_O_type.fasta db_out/2_H_type.fasta \
                       --type nucl \
                       --fasta sample.fasta \
                       --outdir sample_out/ \
                       --threads 2
```

For users own reference sequences files, **seq_typing** requires the construction of the reference database. **seq_typing** will construct the reference DB while analysing the sample's sequences. If many samples will be analysed using the same reference sequences file, a preliminary _seq_typing.py blast_ step is advisable to be run.  

Run **seq_typing** without previous construction of reference database:
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

seq_typing.py assembly --blast references/O_type.fasta references/H_type.fasta \
                       --type nucl \
                       --fasta sample.fasta \
                       --outdir sample_out/ \
                       --threads 2
```

Run **seq_typing** with a preliminary step for reference DB construction (useful when running multiple samples with the same reference sequences file):
```bash
# Preliminary step for reference DB construction.
seq_typing.py blast --blast references/O_type.fasta references/H_type.fasta \
                    --type nucl \
                    --outdir db_out/
# Run seq_typing using created database
seq_typing.py assembly --blast db_out/O_type.fasta db_out/H_type.fasta \
                       --type nucl \
                       --fasta sample.fasta \
                       --outdir sample_out/ \
                       --threads 2
```

### E. coli stx subtyping

A specific script was created for _E. coli_ _stx_ subtyping (**ecoli_stx_subtyping.py**) in order to accommodate the possible existence of _stx2_ paralogs.  
It works very similar to **seq_typing.py**.

#### General usage
```
usage: ecoli_stx_subtyping.py [-h] [--version] {reads,assembly,blast} ...

Gets E. coli stx subtypes

optional arguments:
  -h, --help            Show this help message and exit
  --version             Version information

Subcommands:
  Valid subcommands

  {reads,assembly}
                        Additional help
    reads               reads --help
    assembly            assembly --help
```

**READS**  
Run _ecoli_stx_subtyping.py_ using fastq files.
```
usage: ecoli_stx_subtyping.py reads [-h]
                                    -f /path/to/input/file.fq.gz ...
                                    -r /path/to/reference_sequence.fasta ... | --org stx subtyping
                                    [--stx2covered N] [--stx2identity N]
                                    [-o /path/to/output/directory/] [-j N]
                                    [--mapRefTogether] [--typeSeparator _]
                                    [--extraSeq N] [--minCovPresence N]
                                    [--minCovCall N] [--minGeneCoverage N]
                                    [--minDepthCoverage N] [--minGeneIdentity N]
                                    [--doNotRemoveConsensus] [--debug] [--resume]
                                    [--notClean]

Run ecoli_stx_subtyping.py using fastq files

optional arguments:
  -h, --help            show this help message and exit

Required options:
  -f --fastq /path/to/input/file.fq.gz ...
                        Path to single OR paired-end fastq files. If two files
                        are passed, they will be assumed as being the paired
                        fastq files

Required one of the following options:
  -r --reference 1_virulence_db.stx1_subtyping.fasta 2_virulence_db.stx2_subtyping.fasta
                        Path to stx subtyping reference sequences (if not want to use
                        the ones provided together with seq_typing.py)
  --org stx subtyping   To use stx subtyping reference sequences provided
                        together with seq_typing.py

ecoli_stx_subtyping specific facultative options:
  --stx2covered N       Minimal percentage of sequence covered to consider
                        extra stx2 subtypes (value between [0, 100]) (default: 100)
  --stx2identity N      Minimal sequence identity to consider extra stx2
                        subtypes (value between [0, 100]) (default: 99.5)

General facultative options:
  -o --outdir /path/to/output/directory/
                        Path to the directory where the information will be
                        stored (default: ./)
  -j --threads N        Number of threads to use (default: 1)
  --mapRefTogether      Map the reads against all references together
  --typeSeparator _     Last single character separating the general sequence
                        header from the last part containing the type (default: _)
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
  --minDepthCoverage N  Minimum depth of coverage of target reference sequence
                        to consider a sequence to be present (default: 2)
  --minGeneIdentity N   Minimum percentage of identity of reference sequence
                        covered to consider a gene to be present (value
                        between [0, 100]). One INDEL will be considered as one
                        difference (default: 80)
  --doNotRemoveConsensus
                        Do not remove ReMatCh consensus sequences
  --debug               Debug mode: do not remove temporary files
  --resume              Resume ecoli_stx_subtyping.py reads
```

**ASSEMBLY**  
Run _ecoli_stx_subtyping_ using a fasta file.
```
usage: ecoli_stx_subtyping.py assembly [-h]
                                       -f /path/to/query/assembly_file.fasta
                                       -b /path/to/Blast/db.sequences.file ... -t nucl | --org stx subtyping
                                       [--stx2covered N] [--stx2identity N]
                                       [-o /path/to/output/directory/] [-j N]
                                       [--typeSeparator _] [--minGeneCoverage N]
                                       [--minGeneIdentity N] [--debug]

Run ecoli_stx_subtyping.py using a fasta file. If running multiple samples using the
same DB sequence file, consider use first "seq_typing.py blast"
module.

optional arguments:
  -h, --help            show this help message and exit

Required options:
  -f /path/to/query/assembly_file.fasta, --fasta /path/to/query/assembly_file.fasta
                        Path to fasta file containing the query sequences from
                        which the stx subtypes should be assessed

Required one of the following options:
  -b --blast 1_virulence_db.stx1_subtyping.fasta 2_virulence_db.stx2_subtyping.fasta
                        Path to stx subtyping DB sequence file (if not want to use
                        the ones provided together with seq_typing.py).
                        If Blast DB was already produced (using "seq_typing.py blast"
                        module) only provide the file that do not end with ".n*"
                        something (do not use for example
                        /blast_db.sequences.fasta.nhr). If no Blast DB is
                        found for the DB sequence file, one will be created in
                        --outdir. If more than one Blast DB file is passed, a
                        type for each file will be determined. Give the files
                        in the same order that the type must be determined.
  --org stx subtyping   To use stx subtyping reference sequences provided
                        together with seq_typing.py

Required option for --blast:
  -t --type nucl        Blast DB type (available options: nucl, prot)

ecoli_stx_subtyping specific facultative options:
  --stx2covered 95      Minimal percentage of sequence covered to consider
                        extra stx2 subtypes (value between [0, 100]) (default: 100)
  --stx2identity 95     Minimal sequence identity to consider extra stx2
                        subtypes (value between [0, 100]) (default: 99.5)

General facultative options:
  -o --outdir /path/to/output/directory/
                        Path to the directory where the information will be
                        stored (default: ./)
  -j --threads N        Number of threads to use (default: 1)
  --typeSeparator _     Last single character separating the general sequence
                        header from the last part containing the type (default: _)
  --minGeneCoverage N   Minimum percentage of target reference sequence
                        covered to consider a sequence to be present (value
                        between [0, 100]) (default: 60)
  --minGeneIdentity N   Minimum percentage of identity of reference sequence
                        covered to consider a gene to be present (value
                        between [0, 100]) (default: 80)
  --debug               Debug mode: do not remove temporary files
```

**BLAST**  
To construct stx subtypes Blast DB, proceed as described [here](#assemblies):  
`seq_typing.py blast --org stx subtyping`.

#### Update stx references

An updated stx subtyping reference sequences can be obtained from [VirulenceFinder DB Bitbucket account](https://bitbucket.org/genomicepidemiology/virulencefinder_db). A specific script was created to get the most recent stx reference sequences.
```
usage: get_stx_db.py [-h] [--version]
                     [-o /path/to/output/directory/]

Gets STX sequences from virulencefinder_db to produce a STX subtyping DB.

optional arguments:
  -h, --help            show this help message and exit
  --version             Version information

General facultative options:
  -o --outdir /path/to/output/directory/
                        Path to the directory where the sequences will be
                        stored (default: ./)
```

Usage example
```bash
# Activate Conda environment (when using Python via Conda environment)
source activate seq_typing

get_stx_db.py --outdir /path/output/directory/
```

## Outputs

__seq_typing.report.txt__  
Text file with the typing result. If it was not possible to determine a type for a given reference file, `NT` (for None Typeable) will be returned for that file.

Example of _E. coli_ serotyping (two reference files):  
`O157:H7`  
Example of Dengue virus serotyping and genotyping (only one reference file):  
`3-III`

__seq_typing.report_types.tab__  
Tabular file with detailed results:
* General fields
  * _sequence_type_: type of the results reported. Three values can be found here. `selected` for the reference sequences selected for the reported typing result. `other_probable_type` for other reference sequences that could have been selected because fulfill selection thresholds. `most_likely` for the most likely reference sequences when no reference sequences fulfill selection thresholds.
  * _reference_file_: the reference file where the sequences came from.
  * _type_: the type associated to the reference sequence
  * _sequence_: reference sequences name
  * _sequenced_covered_: percentage of reference sequences covered
  * _coverage_depth_: mean reference sequences depth of coverage of the positions present (1 if assembly was used)
  * _sequence_identity_: percentage identity of reference sequences covered
* Assembly fields (filled with `NA` if reads were used)
  * _query_: name of the provided sequence that had hit with the given reference sequence
  * _q_start_: hit starting position of the provided sequence
  * _q_end_: hit ending position of the provided sequence
  * _s_start_: hit starting position of the reference sequence
  * _s_end_: hit ending position of the reference sequence
  * _evalue_: hit E-value

Example of _E. coli_ serotyping (two reference files) using reads:  

| #sequence_type      | reference_file | type | sequence              | sequenced_covered | coverage_depth     | sequence_identity | query | q_start | q_end | s_start | s_end | evalue |
|---------------------|----------------|------|-----------------------|-------------------|--------------------|-------------------|-------|---------|-------|---------|-------|--------|
| selected            | O_type.fasta   | O26  | wzy_192_AF529080_O26  | 100.0             | 281.95405669599216 | 100.0             | NA    | NA      | NA    | NA      | NA    | NA     |
| selected            | H_type.fasta   | H11  | fliC_269_AY337465_H11 | 99.4546693933197  | 51.76490747087046  | 99.86291980808772 | NA    | NA      | NA    | NA      | NA    | NA     |
| other_probable_type | O_type.fasta   | O26  | wzx_208_AF529080_O26  | 100.0             | 223.3072050673001  | 100.0             | NA    | NA      | NA    | NA      | NA    | NA     |
| other_probable_type | H_type.fasta   | H11  | fliC_276_AY337472_H11 | 98.84117246080436 | 37.52551724137931  | 99.86206896551724 | NA    | NA      | NA    | NA      | NA    | NA     |

Example of Dengue virus serotyping and genotyping (only one reference file) using assembly:  

| #sequence_type      | reference_file                 | type  | sequence                                                 | sequenced_covered | coverage_depth | sequence_identity | query                               | q_start | q_end | s_start | s_end | evalue |
|---------------------|--------------------------------|-------|----------------------------------------------------------|-------------------|----------------|-------------------|-------------------------------------|---------|-------|---------|-------|--------|
| selected            | 1_GenotypesDENV_14-05-18.fasta | 3-III | gb:EU529683|...|Subtype:3-III|Host:Human|seqTyping_3-III | 100.0             | 1              | 99.223            | NODE_1_length_10319_cov_2021.782660 | 138     | 10307 | 10170   | 1     | 0.0    |
| other_probable_type | 1_GenotypesDENV_14-05-18.fasta | 1-V   | gb:GQ868570|...|Subtype:1-V|Host:Human|seqTyping_1-V     | 100.0             | 1              | 99.479            | NODE_2_length_10199_cov_229.028848  | 13      | 10188 | 1       | 10176 | 0.0    |
| other_probable_type | 1_GenotypesDENV_14-05-18.fasta | 4-II  | gb:GQ868585|...|Subtype:4-II|Host:Human|seqTyping_4-II   | 100.0             | 1              | 99.38             | NODE_4_length_10182_cov_29.854132   | 13      | 10173 | 1       | 10161 | 0.0    |

__run.*.log__  
Running log file.  

## Contact

Miguel Machado  
<mpmachado@medicina.ulisboa.pt>
