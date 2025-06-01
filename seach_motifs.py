# -*- coding: utf-8 -*-
import re

# Files
entrada = "intersectFLAG-SCRT2_deduplicado.fa"
saida = "todos_os_motifs_results_sem_filtro_deduplicado.txt"

# Motifs and regex
motifs = {
    "CANNTG": re.compile(r"CA[ACGT][ACGT]TG"),
    "CAACAGGT": re.compile(r"CAACAGGT"),
    "CAGCTTG": re.compile(r"CAGCTTG"),
}

# Searching motifs in a sequence
def buscar_motifs(seq_id, seq, out):
    for nome, padrao in motifs.items():
        for match in padrao.finditer(seq):
            start = match.start()
            match_seq = match.group()
            out.write(f"{seq_id}\t{nome}\t{start}\t{match_seq}\n")

# FASTA processing
with open(entrada) as f_in, open(saida, "w") as f_out:
    seq_id = ""
    seq = ""
    for linha in f_in:
        linha = linha.strip()
        if linha.startswith(">"):
            if seq_id and seq:
                buscar_motifs(seq_id, seq, f_out)
            seq_id = linha[1:]  # remove o >
            seq = ""
        else:
            seq += linha.upper()
    if seq_id and seq:
        buscar_motifs(seq_id, seq, f_out)

