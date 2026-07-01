# ==============================================================================
# PIPELINE DE METABARCODING - PASSO 03: ECOLOGIA MICROBIANA E GRÁFICOS
# Objetivo: Análise de diversidade alfa/beta e composição taxonômica
# ==============================================================================

# 1. Carregar pacotes para ecologia e gráficos
library(phyloseq)
library(ggplot2)
library(vegan)

# 2. Carregar os outputs finais gerados pelo DADA2
load("results/dada2_final_output.RData")

# [Opcional] Criar dados de metadados fictícios para as amostras (ex: Controle vs Tratamento)
# Na vida real, você carregaria um arquivo .csv com as características das suas amostras
sam_data_df <- data.frame(
  Location = rep(c("Solo_Americana", "Solo_Campinas"), each = ncol(seqtab.nochim)/2),
  row.names = rownames(seqtab.nochim)
)

# 3. Construir o objeto Phyloseq (Une todas as tabelas em um único elemento)
OTU <- otu_table(seqtab.nochim, taxa_are_rows = FALSE)
TAX <- tax_table(taxa)
META <- sample_data(sam_data_df)

ps <- phyloseq(OTU, TAX, META)

# 4. Gráfico 1: Diversidade Alfa (Riqueza observada e índice de Shannon)
# Mede a complexidade interna de cada amostra
p_alpha <- plot_richness(ps, x = "Location", measures = c("Observed", "Shannon")) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Diversidade Alfa por Localidade")

ggsave("results/figures/plot_01_alpha_diversity.png", plot = p_alpha, width = 6, height = 4)

# 5. Gráfico 2: Composição Taxonômica (Abundância Relativa a nível de Filo)
# Transforma as contagens absolutas em porcentagem (0 a 100%)
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))
p_bar <- plot_bar(ps_rel, fill = "Phylum") +
  geom_bar(aes(color = Phylum, fill = Phylum), stat = "identity", position = "stack") +
  theme_bw() +
  labs(title = "Composição Taxonômica Relativa (Filo)", y = "Abundância Relativa")

ggsave("results/figures/plot_02_taxa_abundance.png", plot = p_bar, width = 8, height = 5)

# 6. Gráfico 3: Diversidade Beta (Ordenação por NMDS usando distância Bray-Curtis)
# Mostra a similaridade ou diferença espacial entre as comunidades microbianas
ps_ord <- ordinate(ps, "NMDS", "bray")
p_beta <- plot_ordination(ps, ps_ord, color = "Location") +
  geom_point(size = 4) +
  theme_bw() +
  labs(title = "Ordenação NMDS (Diversidade Beta - Bray-Curtis)")

ggsave("results/figures/plot_03_beta_diversity.png", plot = p_beta, width = 6, height = 5)

print("Análise ecológica concluída! Os gráficos foram salvos na pasta results/figures/")
