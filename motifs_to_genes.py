import pandas as pd
import subprocess

# === PART 1: read and conver the motif file ===

# If yout file has 4 columns: coord, motif_type, score, seq
df = pd.read_csv("todos_os_motifs_results_sem_filtro_deduplicado.txt", sep="\t", header=None,
                 names=["coord", "motif_type", "score", "seq"])

# Split coord in chr, start and end
df[["chr", "positions"]] = df["coord"].str.split(":", expand=True)
df[["start", "end"]] = df["positions"].str.split("-", expand=True)
df["start"] = df["start"].astype(int)
df["end"] = df["end"].astype(int)

# Create a BED file with the name being the motif sequence
df_bed = df[["chr", "start", "end", "seq"]]
df_bed.to_csv("motifs.bed", sep="\t", header=False, index=False)

# === PART 2: Run bedtools closest ===

# Make sure that genes.bed already exists
output_file = "motifs_genes_associados.bed"

subprocess.run([
    "bedtools", "closest",
    "-a", "motifs_sorted.bed",
    "-b", "genes_chr_fixed.bed",
    "-D", "ref"  # distância até o gene mais próximo
], stdout=open(output_file, "w"))

# === PART 3: Upload the data in pandas for analysis (optional) ===

cols = [
    "chr", "start", "end", "motif_seq",
    "gene_chr", "gene_start", "gene_end", "gene_name", "score", "strand", "distance"
]
df_out = pd.read_csv(output_file, sep="\t", header=None, names=cols)

# Save
df_out.to_csv("motifs_genes_final_deduplicado.csv", index=False)
