# stx subtyping

Escherichia coli stx subtypes DB

---

## Info

stx subtypes sequences database was originally produced on 27/Nov/2018 from [Virulencefinder Database](https://bitbucket.org/genomicepidemiology/virulencefinder_db) (commit _2159310_) using `get_stx_db.py` script found in `seq_typing/modules/` folder.  

stx A subunits sequences (_stx1A_ and _stx2A_) were updated on 03/Apr/2025 (commit _9638945ea72ec748beded45bb9fe48351eee346f_) using `get_stx_db.py` script.  
stx B subunits sequences were not updated because they don't make part of _Virulencefinder Database_ anymore.  
  
For _stx2_ gene, _stx2a_, _stx2c_ and _stx2d_ variants are grouped together as _stx2acd_
due to the fact that all of these subtypes are the most potent ones to cause HUS and
are difficult to separate from each other by the methods in use right now.
Sequences without subtype were excluded. For sequences having IUPAC codes, a random sequence was produced based on the IUPAC codes.

## Contact

Miguel Machado <mpmachado@medicina.ulisboa.pt>  
Jani Halkilahti <jani.halkilahti@thl.fi>
