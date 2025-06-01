# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt

# Path to BED file
bed_path = "intersectFLAG-SCRT2_deduplicado.bed"

# Read BED and criate PeakID  (chr:start-end)
bed_df = pd.read_csv(bed_path, sep='\t', header=None, names=['chr', 'start', 'end'])
bed_df['PeakID'] = bed_df['chr'].astype(str) + ':' + bed_df['start'].astype(str) + '-' + bed_df['end'].astype(str)
all_peaks = set(bed_df['PeakID'])

# read motif file
motif_df = pd.read_csv('todos_os_motifs_results_sem_filtro.txt', sep='\t', header=None, names=['PeakID', 'Motif', 'Position', 'MatchSequence'])

# Filter PEakIDs in the BED file
motif_df = motif_df[motif_df['PeakID'].isin(all_peaks)]

# Filter unique PEakIDs with motif
peaks_with_motif = set(motif_df['PeakID'])

# Count
with_motif = len(peaks_with_motif)
without_motif = len(all_peaks - peaks_with_motif)

# Organize data
data = pd.Series({
    'With motif': with_motif,
    'None': without_motif
})
percentages = (data / data.sum()) * 100

# Mostrar resultado
result_table = pd.DataFrame({
    'Counts': data,
    'Percentage': percentages
})
print(result_table)

# Pie plot
fig, ax = plt.subplots(figsize=(6, 6))
wedges, texts = ax.pie(
    percentages,
    startangle=140,
    labels=None
)

legend_labels = [f"{label} ({perc:.1f}%)" for label, perc in zip(data.index, percentages)]
plt.legend(
    wedges,
    legend_labels,
    loc="center left",
    bbox_to_anchor=(1, 0, 0.5, 1),
    fontsize=12
)

plt.title('SCRT2 motif presence in unique peaks', fontsize=14)
plt.axis('equal')
plt.tight_layout()

# Save graph
plt.savefig('grafico_scrt2_motif_vs_none.png', dpi=300)
plt.savefig('grafico_scrt2_motif_vs_none.pdf')
print("Gráfico salvo com sucesso!")

# Save CSV
result_table.reset_index(inplace=True)
result_table.rename(columns={'index': 'Category'}, inplace=True)
result_table.to_csv('tabela_motif_vs_none.csv', index=False)
print("Tabela salva como 'tabela_motif_vs_none.csv'")

# Save peaks "None"
peaks_without_motif = all_peaks - peaks_with_motif
none_df = pd.DataFrame({'PeakID': list(peaks_without_motif)})
none_df.to_csv('lista_peaks_None.csv', index=False)
print("Lista de PeakIDs sem motivo salva como 'lista_peaks_None.csv'")
