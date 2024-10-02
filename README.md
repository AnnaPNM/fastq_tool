# fastq_tool
Homework 4 in Python course (BI 2024-2025)

Module contains all basic procedures which could be useful in work with DNA/RNA sequences

## `run_dna_rna_tools`
Different converting procedures for DNA/RNA sequences

**Avaliable procedures**:
- transcribe (return transcribed sequences)
- reverse (return reverse sequences)
- complement (return complement sequences)
- reverse_complement (return reverse complement sequences)

Function accepts any number of DNA/RNA sequences as arguments, *last argument must be the name of procedure*
Function accepts both letter registers

`run_dna_rna_tools` returns a **list** of succesfully processed sequences (if resulting list contains only one element, function returns this element as **str** type)

### Example:
```
run_dna_rna_tools ("AUG", "AAT", "GCCATTT", "reverse")
```
Output:
> All possible procedures were done
> 
> ['GUA', 'TAA', 'TTTACCG']



## `filter_fastq`
Checks the quality of FASTQ reads and filter it

**Arguments**:
- seqs (dict) - dictionary of FASTQ sequences (srt) and its quality (str)
  View: seqs = {"name1":["DNA-seq1", "Quality-DNA-seq1"], "name2":["DNA-seq2", "Quality-DNA-seq2"], ...}
- gc_bounds (tuple) - limit bounds for GC% (default: (0,100)); if only one specified, it will be accepted as upper bound
- length_bounds (tuple) - limit bounds for read length (default: (0, 2**32)); if only one specified, it will be accepted as upper bound
- quality_threshold (float) - limit of mean read quality score (phred33) (default: 0)

`filter_fastq` returns modified dictionary of FASTQ sequences satisfying all conditions

### Example:
```
EXAMPLE_FASTQ = {
    '@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA', 'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD'),
    '@SRX079802': ('ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG', 'BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D')}

filter_fastq (EXAMPLE_FASTQ, gc_bounds = (80), length_bounds = (50, 100), quality_threshold = 32)
```
Output:
> {'@SRX079801': ('ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA',
  'FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD')}




