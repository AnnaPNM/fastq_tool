def if_gc_bounds_ok(seq_, gc_bounds):
    A = seq_.count("A") + seq_.count("a")
    T = seq_.count("T") + seq_.count("t")
    G = seq_.count("G") + seq_.count("g")
    C = seq_.count("C") + seq_.count("c")
    gc_content = (G + C) / (A + T + G + C) * 100
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def if_length_bounds_ok(seq_, length_bounds):
    return length_bounds[0] <= len(seq_) <= length_bounds[1]


def if_quality_threshold_ok(seq_q, quality_threshold):
    seq_q_list = list(seq_q)
    score_list = []
    for symb in seq_q_list:
        score_list.append(ord(symb) - 33)
    return sum(score_list) / len(score_list) >= quality_threshold


def read_fastq(input_fastq):
    with open(input_fastq, "r") as fastq_file:
        seqs = {}
        string_id = 1
        for string in fastq_file:
            if string_id % 4 == 1:
                key_i = string[: len(string) - 1]
            elif string_id % 4 == 2:
                atcg_i = string[: len(string) - 1]
            elif string_id % 4 == 3:
                qual_i = string[: len(string) - 1]
            else:
                seqs[f"{key_i}"] = (atcg_i, string[: len(string) - 1], qual_i)
            string_id += 1
    return seqs


def write_new_fastq(output_fastq, result):
    seq_names = list(result.keys())

    new_fastq = open(output_fastq, "w")

    for seq_ in seq_names:
        new_fastq.write(seq_ + "\n")
        new_fastq.write(result[f"{seq_}"][0] + "\n")
        new_fastq.write(result[f"{seq_}"][2] + "\n")
        new_fastq.write(result[f"{seq_}"][1] + "\n")
    new_fastq.close()
