CHANGE ALL THIS BELOW:

AW: Avoided Words
===

<b>Description</b>: AW is an implementation of a non-trivial linear-time and linear-space algorithm to compute all <b>avoided</b> or all <b>overabundant</b> words in a given DNA or proteins sequence. The definitions used for <b>expectation</b> and <b>variance</b> are described and biologically justified in: 

```
V. Brendel, J.S. Beckmann, and E.N. Trifonov: 
Linguistics of nucleotide sequences: morphology and comparison of vocabularies.
Journal of Biomolecular Structure and Dynamics 4(1), 11-21 (1986).
```

<b>Installation</b>: To compile AW, please follow the instructions given in file INSTALL.

```
 Usage: aw <options>
 Standard (Mandatory):
  -a, --alphabet            <str>     `DNA' for nucleotide  sequences or `PROT'
                                      for protein  sequences. 
  -i, --input-file          <str>     (Multi)FASTA input filename.
  -o, --output-file         <str>     Output filename.
  -t, --threshold           <dbl>     The threshold.
 Optional:
  -w, --words-class         <int>     `0' to check for AVOIDED words or `1' to
                                      check for OVERABUNDANT (default: 0).
  -k, --length              <int>     Fixed length of words (default: search all).
  -A, --absent              <int>     `1' to check also for ABSENT AVOIDED words
                                      or `0' otherwise (default: 0).
  -r, --reverse             <int>     `1' to check for the reverse complement or
                                      `0' otherwise (default: 0).
```

<b>Citations</b>:

```
Y. Almirantis, P. Charalampopoulos, J. Gao, C. S. Iliopoulos, M. Mohamed, S. P. Pissis, D. Polychronopoulos: 
Optimal Computation of Avoided Words. 
WABI 2016: 1-13.

Y. Almirantis, P. Charalampopoulos, J. Gao, C. S. Iliopoulos, M. Mohamed, S. P. Pissis, D. Polychronopoulos: 
Optimal Computation of Overabundant Words. 
WABI 2017: 4:1-4:14.
```
<b>License</b>: GNU GPLv3 License; Copyright (C) 2016 Jia Gao and Solon P. Pissis
