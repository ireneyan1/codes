# Pacotes necessários
library(ChIPseeker)
library(GenomicRanges)
library(rtracklayer)
library(GenomicFeatures)
library(GenomeInfoDb)
library(org.Gg.eg.db)
library(dplyr)
library(ggplot2)
library(seaborn)  # para compatibilidade com sns (caso esteja usando reticulate com Python)

# Caminhos de entrada
bed_file <- "/home/geninfo/vbotezelli/intersectFLAG-SCRT2_deduplicado.bed"
gtf_file <- "/home/geninfo/vbotezelli/Gallus_gallus.GRCg6a.106.gtf"

# Criar TxDb a partir do GTF
UCSCgg6_txdb <- makeTxDbFromGFF(gtf_file, format = 'gtf')

# Importar BED e ajustar estilo dos cromossomos
Sctr_Flagcs_ls <- import(bed_file, format = 'bed')
Sctr_Flagcs_gr <- GRanges(Sctr_Flagcs_ls)
seqlevelsStyle(Sctr_Flagcs_gr) <- seqlevelsStyle(UCSCgg6_txdb)

# Anotar os picos
Sctr_Flagcs_gr.annot <- annotatePeak(Sctr_Flagcs_gr, TxDb = UCSCgg6_txdb, annoDb = 'org.Gg.eg.db')
annot_df <- as.data.frame(Sctr_Flagcs_gr.annot)
write.csv(annot_df, 'annotated_peaks.csv', row.names = FALSE)

# Agrupar categorias: Distal Intergenic, Promoter, Exon, Intron, Other
annot_df$annotation_grouped <- dplyr::case_when(
  grepl('Distal Intergenic', annot_df$annotation, ignore.case = TRUE) ~ 'Distal Intergenic',
  grepl('Promoter', annot_df$annotation, ignore.case = TRUE) ~ 'Promoter',
  grepl('Exon', annot_df$annotation, ignore.case = TRUE) ~ 'Exon',
  grepl('Intron', annot_df$annotation, ignore.case = TRUE) ~ 'Intron',
  TRUE ~ 'Other'
)

# Contagem por grupo
grouped_counts <- annot_df %>%
  dplyr::count(annotation_grouped, name = 'Contagem') %>%
  dplyr::mutate(
    Percentual = round(100 * Contagem / sum(Contagem), 1),
    Legenda = paste0(annotation_grouped, ' (', Percentual, '%)')
  )

# Gráfico de pizza com legenda ao lado
png('custom_plot_annotation_pie.png', width = 800, height = 600, res = 150)
pie(grouped_counts$Contagem, labels = NA, main = 'Peak Distribution by Category')
legend('right', legend = grouped_counts$Legenda, bty = 'n', xpd = TRUE)
dev.off()

# Gráfico de barras horizontais
grouped_counts_sorted <- grouped_counts[order(grouped_counts$Contagem), ]
png('custom_plot_annotation_bar.png', width = 800, height = 600, res = 150)
barplot(grouped_counts_sorted$Contagem, names.arg = grouped_counts_sorted$annotation_grouped, horiz = TRUE, las = 1, col = 'skyblue', xlab = 'Peak Count')
dev.off()

# Gráfico empilhado
png('custom_plot_annotation_stacked.png', width = 800, height = 200, res = 150)
barplot(height = grouped_counts$Contagem, names.arg = NA, col = rainbow(nrow(grouped_counts)), horiz = TRUE, axes = FALSE)
legend('topright', legend = grouped_counts$Legenda, fill = rainbow(nrow(grouped_counts)), bty = 'n')
dev.off()