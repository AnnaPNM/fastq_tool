def is_dna(seq_):
    return set(list(seq_)) <= set(["A", "T", "C", "G", "a", "t", "c", "g"])


def is_rna(seq_):
    return set(list(seq_)) <= set(["A", "U", "C", "G", "a", "u", "c", "g"])


def transcribe_seq(seq_):
    res = []
    list_seq = list(seq_)
    for n in list_seq:
        if n == "T":
            res.append("U")
        elif n == "t":
            res.append("u")
        else:
            res.append(n)
    return "".join(res)


def reverse_seq(seq_):
    res = []
    list_seq = list(seq_)
    for i in range(len(list_seq) - 1, -1, -1):
        res.append(list_seq[i])
    return "".join(res)


def complement_seq(seq_):
    res = []
    list_seq = list(seq_)
    for n in list_seq:
        if n == "T":
            res.append("A")
        elif n == "t":
            res.append("a")
        elif n == "A":
            if "U" in list_seq or "u" in list_seq:
                res.append("U")
            else:
                res.append("T")
        elif n == "a":
            if "U" in list_seq or "u" in list_seq:
                res.append("u")
            else:
                res.append("t")
        elif n == "C":
            res.append("G")
        elif n == "c":
            res.append("g")
        elif n == "G":
            res.append("C")
        elif n == "g":
            res.append("c")
        elif n == "U":
            res.append("A")
        elif n == "u":
            res.append("a")

    return "".join(res)


def reverse_complement_seq(seq_):
    return reverse_seq(complement_seq(seq_))
