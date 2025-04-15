from abc import ABC, abstractmethod
from typing import Union

import Bio
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord


class BiologicalSequence(ABC):
    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def is_correct_seq(self):
        pass


class NucleicAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, val):
        return self.sequence[val]

    def __repr__(self):
        return f"NucleicAcidSequence {self.sequence}"

    def __str__(self):
        return str(self.sequence)

    def is_correct_seq(self):
        return set(self.sequence) <= {"A", "T", "C", "G", "U", "a", "t", "c", "g", "u"}

    def complement(self):
        res = []
        list_seq = list(self.sequence)

        if "U" in list_seq or "u" in list_seq:
            complement_list = {
                "A": "U",
                "U": "A",
                "C": "G",
                "G": "C",
                "a": "u",
                "u": "a",
                "c": "g",
                "g": "c",
            }
        else:
            complement_list = {
                "A": "T",
                "T": "A",
                "C": "G",
                "G": "C",
                "a": "t",
                "t": "a",
                "c": "g",
                "g": "c",
            }

        for nucl in list_seq:
            res.append(complement_list[nucl])

        return self.__class__("".join(res))

    def reverse(self):
        res = []
        list_seq = list(self.sequence)
        for i in range(len(list_seq) - 1, -1, -1):
            res.append(list_seq[i])
        return self.__class__("".join(res))

    def reverse_complement(self):
        return self.__class__((self.complement()).reverse())


class RNASequence(NucleicAcidSequence):
    pass


class DNASequence(NucleicAcidSequence):
    def transcribe(self):
        res = []
        list_seq = list(self.sequence)

        DNA_to_RNA = {
            "A": "U",
            "T": "A",
            "C": "G",
            "G": "C",
            "a": "u",
            "t": "a",
            "c": "g",
            "g": "c",
        }

        for nucl in list_seq:
            res.append(DNA_to_RNA[nucl])
        return RNASequence("".join(res))


class AminoAcidSequence(BiologicalSequence):
    def __init__(self, sequence: str):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, val):
        return self.sequence[val]

    def __repr__(self):
        return f"AminoAcidSequence {self.sequence}"

    def __str__(self):
        return str(self.sequence)

    def is_correct_seq(self):
        return set(self.sequence) <= {
            "A",
            "R",
            "N",
            "D",
            "C",
            "E",
            "Q",
            "G",
            "H",
            "I",
            "L",
            "K",
            "M",
            "F",
            "P",
            "S",
            "T",
            "W",
            "Y",
            "V",
            "a",
            "r",
            "n",
            "d",
            "c",
            "e",
            "q",
            "g",
            "h",
            "i",
            "l",
            "k",
            "m",
            "f",
            "p",
            "s",
            "t",
            "w",
            "y",
            "v",
        }

    def get_mass_dalton(self):
        """Method get_mass_dalton()
        Returns molecular mass of given amino-acid sequence in Daltons
        """

        mass = 0
        amino_acid_weights = {
            "A": 89.09,  # Alanine
            "R": 174.20,  # Arginine
            "N": 132.12,  # Asparagine
            "D": 133.10,  # Aspartic acid
            "C": 121.15,  # Cysteine
            "E": 147.13,  # Glutamic acid
            "Q": 146.14,  # Glutamine
            "G": 75.07,  # Glycine
            "H": 155.16,  # Histidine
            "I": 131.17,  # Isoleucine
            "L": 131.17,  # Leucine
            "K": 146.19,  # Lysine
            "M": 149.21,  # Methionine
            "F": 165.19,  # Phenylalanine
            "P": 115.13,  # Proline
            "S": 105.09,  # Serine
            "T": 119.12,  # Threonine
            "W": 204.23,  # Tryptophan
            "Y": 181.19,  # Tyrosine
            "V": 117.15,
        }  # Valine
        for acid in list(self.sequence):
            mass += amino_acid_weights[acid]

        return mass


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

    seqs = SeqIO.parse(input_fastq, "fastq")

    result = []

    if not (isinstance(gc_bounds, tuple)):
        gc_bounds = (0, gc_bounds)
    if not (isinstance(length_bounds, tuple)):
        length_bounds = (0, length_bounds)

    for seq_ in seqs:
        mean_quality = sum(seq_.letter_annotations["phred_quality"]) / len(
            seq_.letter_annotations["phred_quality"]
        )
        gc_content = SeqUtils.gc_fraction(seq_.seq) * 100
        seq_length = len(seq_.seq)

        if (
            mean_quality >= quality_threshold
            and gc_bounds[0] <= gc_content <= gc_bounds[1]
            and length_bounds[0] <= seq_length <= length_bounds[1]
        ):
            result.append(SeqRecord(seq_.seq, id=seq_.id, description=seq_.description))

    SeqIO.write(result, output_fastq, "fasta")
