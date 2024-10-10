def convert_multiline_fasta_to_oneline(input_fasta: str,
                                       output_fasta: str = None):
    if output_fasta is None:
        output_fasta = input_fasta.split(".")
        output_fasta = "".join(output_fasta[: len(output_fasta) - 1])
        output_fasta = output_fasta + "_oneline.fasta"

    with open(input_fasta, "r") as fasta_file:
        new_fasta = open(output_fasta, "w")
        new_line = ""
        for line in fasta_file:
            if line[0] == ">":
                if new_line != "":
                    new_fasta.write(new_line + "\n")
                new_fasta.write(line)
                new_line = ""
            else:
                new_line = new_line + line[: len(line) - 1]
        new_fasta.write(new_line)
        new_fasta.close()


def parse_blast_output(input_file: str, output_file: str = None):
    if output_file is None:
        output_file = input_file.split(".")
        output_file = "".join(output_file[: len(output_file) - 1])
        output_file = output_file + "_proteins.txt"

    with open(input_file, "r") as file:
        proteins = []

        line_count = 0

        for line in file:
            if line == "Sequences producing significant alignments:\n":
                line_count = 1
            elif line_count == 1:
                line_count = 2
            elif line_count == 2:
                line_count = 3
            elif line_count == 3:
                proteins.append(line.split("  ")[0])
                line_count = 0

        proteins.sort()
        proteins = "\n".join(proteins)

        new_file = open(output_file, "w")
        new_file.write(proteins)
        new_file.close()


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: tuple,
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = None,
):
    if output_fasta is None:
        output_fasta = input_gbk.split(".")
        output_fasta = "".join(output_fasta[: len(output_fasta) - 1])
        output_fasta = output_fasta + "_output.fasta"

    with open(input_gbk, "r") as gbk_file:
        genes_anno = {}

        flag = False
        new_seq = "no_seq"

        gene_number = 1

        for line in gbk_file:
            if "/gene=" in line:
                new_gene = line.split('"')[1]
                flag = True

            if flag and not ("/gene=" in line):
                if "/translation=" in line:
                    new_seq = line.split('="')[1]
                    new_seq = new_seq[: len(new_seq) - 1]
                if (
                    new_seq[-1] != '"'
                    and new_seq != "no_seq"
                    and not ("/translation=" in line)
                ):
                    seq_add = line.split("  ")[-1]
                    new_seq = new_seq + seq_add[1 : len(seq_add) - 1]
                elif new_seq[-1] == '"':
                    new_seq = new_seq[: len(new_seq) - 1]
                    genes_anno[f"{new_gene}"] = (new_seq, gene_number)
                    gene_number = gene_number + 1
                    flag = False
                    new_seq = "no_seq"

    genes_anno_keys = list(genes_anno.keys())

    new_fasta = open(output_fasta, "w")

    for gene in genes:
        match_genes = [s for s in genes_anno_keys if gene in s]
        if match_genes == []:
            print(f"Error! Gene {gene} is not found")
        else:
            for match_gene in match_genes:
                gene_number = genes_anno[f"{match_gene}"][1]
                if gene_number - n_before < 0:
                    print(f"Error! There is no {n_before} genes before {gene}")
                elif gene_number + n_after > len(genes_anno):
                    print(f"Error! There is no {n_after} genes after {gene}")
                else:
                    n_start = gene_number - n_before
                    n_end = gene_number + n_after
                    for i in range(n_start, n_end + 1):
                        gene_i = genes_anno_keys[i - 1]
                        new_fasta.write(
                            f">{gene_i} ({match_gene} range:[-{n_before};+{n_after}])\n"
                        )
                        new_fasta.write(genes_anno[f"{gene_i}"][0] + "\n")

    new_fasta.close()
