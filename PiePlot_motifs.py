# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt

# Read motif file
df = pd.read_csv('todos_os_motifs_results_sem_filtro_deduplicado.txt', sep='\t', header=None, names=['PeakID', 'Motif', 'Position', 'MatchSequence'])

# Map names
motif_name_map = {
    'CANNTG': 'Ebox',
    'CAGCTTG': 'vCESbox',
    'CAACAGGT': 'iCESbox'
}
df['Motif'] = df['Motif'].map(motif_name_map)

# Count and percentage
motif_counts = df['Motif'].value_counts()
motif_percentages = (motif_counts / motif_counts.sum()) * 100

# Total and table
total_matches = motif_counts.sum()
print("Total de matches encontrados: {}".format(total_matches))

result_table = pd.DataFrame({
    'Counts': motif_counts,
    'Percentage': motif_percentages
})
print(result_table)

# Pie plot
fig, ax = plt.subplots(figsize=(8,8))

wedges, texts = ax.pie(
    motif_percentages,
    startangle=140,
    labels=None
)

legend_labels = [f"{name} ({perc:.1f}%)" for name, perc in zip(motif_percentages.index, motif_percentages)]
plt.legend(
    wedges,
    legend_labels,
    loc="center left",
    bbox_to_anchor=(1, 0, 0.5, 1),
    fontsize=12,
    title=None
)

plt.title('SCRT2 binding sites distribution', fontsize=16)
plt.axis('equal')
plt.tight_layout()

# Save as PNG and PDF
plt.savefig('grafico_scrt2_pizza_legenda.png', dpi=300)
plt.savefig('grafico_scrt2_pizza_legenda.pdf')
