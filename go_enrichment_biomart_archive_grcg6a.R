# Instale os pacotes necessários, se não tiver
# BiocManager::install(c('biomaRt', 'clusterProfiler', 'org.Gg.eg.db', 'enrichplot'))
library(biomaRt)
library(clusterProfiler)
library(org.Gg.eg.db)
library(enrichplot)
library(ggplot2)

# Carregar os genes da anotação
annot_df <- read.csv("annotated_peaks.csv")
genes_ensembl <- unique(na.omit(annot_df$geneId))

# Conectar ao Ensembl antigo (abril 2020), compatível com GRCg6a
ensembl_old <- useMart(
  'ENSEMBL_MART_ENSEMBL',
  dataset = 'ggallus_gene_ensembl',
  host = 'https://apr2020.archive.ensembl.org'
)

# Obter correspondência ENSEMBL → ENTREZ
conversion <- getBM(
  attributes = c('ensembl_gene_id', 'entrezgene_id'),
  filters = 'ensembl_gene_id',
  values = genes_ensembl,
  mart = ensembl_old
)

# Limpar e extrair ENTREZ IDs válidos
gene_entrez <- na.omit(unique(conversion$entrezgene_id))

# Rodar enriquecimento GO
ego <- enrichGO(
  gene = gene_entrez,
  OrgDb = org.Gg.eg.db,
  keyType = 'ENTREZID',
  ont = 'BP',
  pAdjustMethod = 'BH',
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.2,
  readable = TRUE
)

# Exportar resultados
ego_df <- as.data.frame(ego)
write.csv(ego_df, 'go_enrichment_results_biomart_archive.csv', row.names = FALSE)

# Gráficos
png('go_enrichment_barplot_biomart_archive.png', width = 1000, height = 800, res = 150)
barplot(ego, showCategory = 15, title = 'GO Enrichment - BP (GRCg6a, BioMart 2020)')
dev.off()

pdf('go_enrichment_barplot_biomart_archive.pdf', width = 10, height = 8)
barplot(ego, showCategory = 15, title = 'GO Enrichment - BP (GRCg6a, BioMart 2020)')
dev.off()

png('go_enrichment_dotplot_biomart_archive.png', width = 1000, height = 800, res = 150)
dotplot(ego, showCategory = 15, title = 'GO Enrichment - BP (GRCg6a, BioMart 2020)')
dev.off()

pdf('go_enrichment_dotplot_biomart_archive.pdf', width = 10, height = 8)
dotplot(ego, showCategory = 15, title = 'GO Enrichment - BP (GRCg6a, BioMart 2020)')
dev.off()

library(ggplot2)
library(viridis)
library(stringr)
# Extrair e preparar os dados
df <- as.data.frame(ego)[1:6, ]
df$Description <- str_to_sentence(df$Description)  # Letra maiúscula no início
df$Description <- factor(df$Description, levels = rev(df$Description))  # manter ordem no eixo

# Definir cores da paleta inferno
inferno_colors <- viridis::inferno(n = 6)

# Criar o gráfico
ggplot(df, aes(x = Description, y = Count, fill = Description)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = inferno_colors) +
  labs(title = "GO Enrichment - BP (GRCg6a, BioMart 2020)", x = NULL, y = "Count") +
  theme_minimal(base_size = 16) +  # base_size controla a maioria dos textos
  theme(
    legend.position = "none",
    axis.text.y = element_text(size = 16, color = "black", face = "plain"),  # fonte preta, maior
    axis.text.x = element_text(size = 14, color = "black"),
    axis.title.x = element_text(size = 14, color = "black"),
    plot.title = element_text(size = 18, face = "bold", color = "black"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
ggsave("go_enrichment_barplot_biomart_archive.png", width = 12, height = 8, dpi = 150)
