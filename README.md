# SequenceAlignment
Provides a command, `gsa`, for global alignment of two sequences (in FASTA format) using a scoring matrix.

Two example sequences are provided:
* TestData/A2ASS6.fasta: TITIN_MOUSE
* TestData/Q8WZ42.fasta: TITIN_HUMAN

An example scoring matrix is also provided:
* TestData/blosum62.txt

### Usage

```
gsa [program options] <mode [options]>

[program options]:
    -h, --help          This help
    -b, --bottom FILE   file with the bottom sequence (default:
                        TestData/Q8WZ42.fasta)
    -t, --top FILE      file with the top sequence (default:
                        TestData/A2ASS6.fasta)
    -s, --scores FILE   The file with the score matrix (default:
                        TestData/blosum62.txt)
    -g, --gip INT       The gap insertion penalty (default: -10)
    -a, --align FILE    The file that will store the alignment

<mode [options]>:
    info                      Shows information about the matrix
                              score and the two sequences
    bl                        Baseline (only score, no alignment)
    qs                        Quadratic-Space DP algorithm
    hb <rec_threshold>        Hirschberg algorithm
    kcol <k> <rec_threshold>  K-Col algorithm
```

### Usage examples
```
gsa bl
gsa qs
gsa hb 20000
gsa kcol 32 20000
```