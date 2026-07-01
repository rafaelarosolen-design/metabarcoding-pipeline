# ==============================================================================
# PIPELINE DE METABARCODING - PASSO 02: TAXONOMIA E TABELA DE ASVs
# Objetivo: Merge de reads, remoção de quimeras e atribuição taxonômica
# ==============================================================================

library(dada2)

# 1. Carregar os dados processados no Passo 01
load("results/dada2_partial_output.RData")

# 2. Unir as leituras Forward e Reverse (Merge)
# Alinha as duas fitas pareadas para formar a sequência completa do amplicon
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

# 3. Construir a tabela de sequências (Sequence Table)
# É a matriz equivalente à "OTU table", mostrando a abundância de cada ASV por amostra
seqtab <- makeSequenceTable(mergers)

# 4. Remover Quimeras
# Quimeras são artefatos de PCR onde duas sequências distintas se unem erroneamente.
# O DADA2 compara as sequências e remove esses erros biológicos.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

# 5. Atribuição Taxonômica (Classificação)
# NOTA: Você precisa baixar o banco de dados de referência (ex: SILVA para 16S, UNITE para ITS)
# Substitua "banco_de_referencia.fasta.gz" pelo caminho real do arquivo do banco.
taxa <- assignTaxonomy(seqtab.nochim, "data/reference_db/banco_de_referencia.fasta.gz", multithread=TRUE)

# 6. Salvar os resultados finais estruturados
# Salva em formato RData para abrir direto no RStudio
save(seqtab.nochim, taxa, file = "results/dada2_final_output.RData")

# Exportar também em arquivos de texto (.csv) para usar em outros programas se quiser
write.csv(as.data.frame(t(seqtab.nochim)), "results/asv_table.csv")
write.csv(as.data.frame(taxa), "results/taxonomy_table.csv")

print("Pipeline do DADA2 finalizado! Tabelas prontas para análise estatística.")
