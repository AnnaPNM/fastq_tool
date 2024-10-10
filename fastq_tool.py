from modules.dna_rna_tools_extended import (
    is_dna,
    is_rna,
    transcribe_seq,
    reverse_seq,
    complement_seq,
    reverse_complement_seq,
)
from modules.filter_fastq_extended import (
    if_gc_bounds_ok,
    if_length_bounds_ok,
    if_quality_threshold_ok,
    read_fastq,
    write_new_fastq,
)

from typing import Union


def run_dna_rna_tools(*args: str) -> Union[str, list]:
    """Function run_dna_rna_tools
    Args: DNA/RNA sequences (str),
    last argument - name of procedure of interest (str)
        Possible procedures:
            transcribe - return transcribed sequences
            reverse - return reverse sequences
            complement - return complement sequences
            reverse_complement - return reverse complement sequences
    Returns:
        If all the sequences were succesfully processed,
        run_dna_rna_tools returns a list of them
        If resulting list contains only one element,
        run_dna_rna_tools returns this element as str type
    """
    procedure = args[-1]
    args_seq = args[: len(args) - 1]
    result = []

    if procedure == "transcribe":
        for seq_ in args_seq:
            if is_dna(seq_) or is_rna(seq_):
                result.append(transcribe_seq(seq_))

    elif procedure == "reverse":
        for seq_ in args_seq:
            if is_dna(seq_) or is_rna(seq_):
                result.append(reverse_seq(seq_))

    elif procedure == "complement":
        for seq_ in args_seq:
            if is_dna(seq_) or is_rna(seq_):
                result.append(complement_seq(seq_))

    elif procedure == "reverse_complement":
        for seq_ in args_seq:
            if is_dna(seq_) or is_rna(seq_):
                result.append(reverse_complement_seq(seq_))
    else:
        print("Unknown procedure:", procedure, "!!!")

    if result == []:
        print("All arguments are not DNA/RNA")
    else:
        print("All possible procedures were done")
        if len(result) == 1:
            return result[0]
        else:
            return result


def filter_fastq(
    input_fastq: str,
    output_fastq: str,
    gc_bounds: Union[tuple, int, float] = (0, 100),
    length_bounds: Union[tuple, int] = (0, 2**32),
    quality_threshold: Union[float, int] = 0,
):
    """Function filter_fastq
    Args:
        input_fastq (str) - path to fastq file
        output_fastq (str) - path (and name) to resulting fastq file
        gc_bounds (tuple) - limit bounds for GC% (default: (0,100));
        if only one specified, it will be accepted as upper bound
        length_bounds (tuple) - limit bounds for read length
            (default: (0, 2**32));
            if only one specified, it will be accepted as upper bound
        quality_threshold (float) - limit of
            mean read quality score (phred33) (default: 0)
    Returns: None
    Create new file with filtered reads (output_fastq)
    """

    seqs = read_fastq(input_fastq)

    result = dict()

    seq_names = list(seqs.keys())

    if type(gc_bounds) != tuple:
        gc_bounds = (0, gc_bounds)
    if type(length_bounds) != tuple:
        length_bounds = (0, length_bounds)

    for seq_name_i in seq_names:
        seq_ = seqs[f"{seq_name_i}"][0]
        seq_q = seqs[f"{seq_name_i}"][1]
        seq_q_id = seqs[f"{seq_name_i}"][2]
        if (
            if_gc_bounds_ok(seq_, gc_bounds)
            and if_length_bounds_ok(seq_, length_bounds)
            and if_quality_threshold_ok(seq_q, quality_threshold)
        ):
            result[f"{seq_name_i}"] = (seq_, seq_q, seq_q_id)

    write_new_fastq(output_fastq, result)
