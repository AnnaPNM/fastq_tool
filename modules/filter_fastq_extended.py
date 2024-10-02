def if_gc_bounds_OK(seq_, gc_bounds):
    A = seq_.count("A") + seq_.count("a")
    T = seq_.count("T") + seq_.count("t")
    G = seq_.count("G") + seq_.count("g")
    C = seq_.count("C") + seq_.count("c")
    gc_content = (G + C) / (A + T + G + C) * 100
    return gc_bounds[0] <= gc_content <= gc_bounds[1]


def if_length_bounds_OK(seq_, length_bounds):
    return length_bounds[0] <= len(seq_) <= length_bounds[1]


def if_quality_threshold_OK(seq_q, quality_threshold):
    seq_q_list = list(seq_q)
    score_list = []
    for symb in seq_q_list:
        score_list.append(ord(symb) - 33)
    return sum(score_list) / len(score_list) >= quality_threshold
