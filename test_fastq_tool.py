import os
# import pytest

# from fastq_tool import (
#     BiologicalSequence,
#     NucleicAcidSequence,
#     RNASequence,
#     DNASequence,
#     AminoAcidSequence,
# )

from fastq_tool import filter_fastq


def test_filtering():
    filter_fastq(
        input_fastq="example_fastq_test.fastq",
        output_fastq="output_example.fastq",
        gc_bounds=30,
        length_bounds=50,
        quality_threshold=30,
    )

    number_of_strings = 0

    with open("output_test.fastq", "r") as f:
        for string in f:
            number_of_strings += 1

    assert number_of_strings == 2


def test_correct_writing_format():
    filter_fastq(
        input_fastq="example_fastq_test.fastq",
        output_fastq="output_example.fastq",
        gc_bounds=30,
        length_bounds=50,
        quality_threshold=30,
    )

    test_result = []

    with open("output_test.fastq", "r") as f:
        for string in f:
            test_result.append(string)

    assert (
        test_result[0] == ">SRX079804:1:SRR292678:1:1101:703304:703304 2:N:0:1 BH:ok\n"
    )
    assert test_result[1] == "TAATAATATAAATTGCTTCTGCTTCTAATTTATCAAGATGTGATAA\n"


def test_empty_output_due_to_high_quality_threshold():
    filter_fastq(
        input_fastq="example_fastq_test.fastq",
        output_fastq="output_high_quality.fastq",
        gc_bounds=0,
        length_bounds=2**32,
        quality_threshold=100,
    )
    with open("output_high_quality.fastq", "r") as f:
        content = f.read().strip()
    assert content == ""


def test_write_fastq_exists():
    filter_fastq(
        input_fastq="example_fastq_test.fastq",
        output_fastq="output_example.fastq",
        gc_bounds=30,
        length_bounds=50,
        quality_threshold=30,
    )

    assert os.path.exists("output_example.fastq")


def test_file_not_found_error():
    try:
        filter_fastq(
            input_fastq="not_existing_file.fastq",
            output_fastq="output.fastq",
            gc_bounds=30,
            length_bounds=50,
            quality_threshold=30,
        )
    except Exception as e:
        assert type(e).__name__ == "FileNotFoundError"


def test_gc_bounds_incorrect_type_error():
    try:
        filter_fastq(
            input_fastq="example_fastq_test.fastq",
            output_fastq="output.fastq",
            gc_bounds="NOT_INT",
            length_bounds=50,
            quality_threshold=30,
        )
    except Exception as e:
        assert type(e).__name__ == "TypeError"


def test_length_bounds_incorrect_type_error():
    try:
        filter_fastq(
            input_fastq="example_fastq_test.fastq",
            output_fastq="output.fastq",
            gc_bounds=30,
            length_bounds="NOT_INT",
            quality_threshold=30,
        )
    except Exception as e:
        assert type(e).__name__ == "TypeError"


def test_quality_threshold_incorrect_type_error():
    try:
        filter_fastq(
            input_fastq="example_fastq_test.fastq",
            output_fastq="output.fastq",
            gc_bounds=30,
            length_bounds=50,
            quality_threshold="NOT_INT",
        )
    except Exception as e:
        assert type(e).__name__ == "TypeError"
