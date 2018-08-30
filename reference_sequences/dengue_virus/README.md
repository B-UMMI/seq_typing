# Dengue virus typing

Dengue virus types DB (serotypes and genotypes).

---

## Info

Dengue virus types sequences database obtained on 09/May/2018 from the NIAID Virus Pathogen Database and Analysis Resource ([ViPR](http://www.viprbrc.org/)).  
The selection criteria for the search were as follows:
  * complete sequence of ORF
  * Human host only
  * collection year (1950-2018)

Data available from all countries was included. 
 
Sequences were subsequently cleaned with `check_sequences.py` script. Sequences without explicit _"Subtype"_ were excluded. Those that have weighted IUPAC codes (weighted by the number of nucleotide possibilities coded for each IUPAC code) frequency >= 0.002 were also excluded. For the other sequences having IUPAC codes, a random sequence was produced based on the IUPAC codes. Finally, redundancy was removed by collapsing identical sequences and contained sequences into a single entry.

## Contact

Miguel Machado: <mpmachado@medicina.ulisboa.pt>  
Erley Lizarazo: <e.f.lizarazo.forero@umcg.nl>
