# fastq_tool
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
- input_fastq (str) - path to the input fastq file
- output_fastq (str) - path (and name) of the output fastq file (if None - `filter_fastq` creates an output file in the working directory) (default: None)
- gc_bounds (tuple) - limit bounds for GC% (default: (0,100)); if only one specified, it will be accepted as upper bound
- length_bounds (tuple) - limit bounds for read length (default: (0, 2**32)); if only one specified, it will be accepted as upper bound
- quality_threshold (float) - limit of mean read quality score (phred33) (default: 0)

`filter_fastq` creates a new file containing FASTQ sequences satisfying all conditions

### Usage:
```
filter_fastq (input_fastq = "file.fastq", gc_bounds = (80), length_bounds = (50, 100), quality_threshold = 32)
```


# bio_files_processor
Module contains procedures for work with "bio"-files

## `convert_multiline_fasta_to_oneline`
Create new one-line (1 sequence in the 1 string) FASTA file from multi-line FASTA file

**Arguments**:
- input_fasta (str) - path to the input FASTA file
- output_fasta (str) - path (and name) of the output FASTA file (if None - `convert_multiline_fasta_to_oneline` creates an output file in the working directory) (default: None)

### Usage:
```
convert_multiline_fasta_to_oneline(input_fasta = 'example_multiline_fasta.fasta')
```

## `parse_blast_output`
Create new file with sorted list of proteins from BLAST-result file

**Arguments**:
- input_file (str) - path to the input txt file with BLAST result
- output_file (str) - path (and name) of the output file (if None - `parse_blast_output` creates an output file in the working directory) (default: None)

### Usage:
```
parse_blast_output(input_file = 'example_blast_results.txt')
```

## `select_genes_from_gbk_to_fasta`
Create new FASTA file with protein sequences of genes of interest in the particular interval from *.gbk file

**Arguments**:
- input_gbk (str) - path to the input *.gbk file with BLAST result
- genes (tuple) - names of genes of interest
- n_before (int) - number of genes before gene of interest (default: 1)
- n_after (int) - number of genes after gene of interest (default: 1)
- output_fasta (str) - path (and name) of the output file (if None - `select_genes_from_gbk_to_fasta` creates an output file in the working directory) (default: None)

If there are several lines with information about gene, `select_genes_from_gbk_to_fasta` will find nearest genes for all of them 
For example, your gene of interest is btuD, n_before = 4 and n_after = 3. In *.gbk file there are three btuD entries: btuD_1, btuD_2, btuD_3. In this situation  `select_genes_from_gbk_to_fasta` will find protein sequences of 4 genes before each of (btuD_1, btuD_2, btuD_3) and of 3 genes after each of (btuD_1, btuD_2, btuD_3) including protein sequences of (btuD_1, btuD_2, btuD_3).

### Usage:
```
select_genes_from_gbk_to_fasta(input_gbk = 'example_gbk.gbk', 
                               genes = ('dtpD', 'clpA', 'rhsD'),
                               n_before = 1, n_after = 3)
```
