# ==============================================================================
# PIPELINE DE METABARCODING - PASSO 01: PROCESSAMENTO COM DADA2
# Objetivo: Filtragem, modelagem de erros e inferência de ASVs
# ==============================================================================

# 1. Carregar o pacote necessário
library(dada2)

# 2. Definir caminhos das pastas (Baseado na estrutura do nosso repositório)
path <- "data/raw_data" # Pasta com os arquivos .fastq brutos
path_filtered <- "data/filtered_data" # Pasta onde salvaremos os arquivos limpos

# Listar arquivos de leitura (Reads) Forward (R1) e Reverse (R2)
# Exemplo padrão de nomenclatura: amostra_R1.fastq / amostra_R2.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))

# Extrair nomes das amostras
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

# 3. Definir caminhos para os arquivos filtrados que serão gerados
filtFs <- file.path(path_filtered, paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path_filtered, paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

# 4. Filtragem e Corte (Filter and Trim)
# NOTA: Os parâmetros truncLen devem ser ajustados conforme a qualidade dos seus dados!
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(240,160), # Corta as leituras F e R nesses comprimentos
                     maxN=0,              # DADA2 não aceita bases "N"
                     maxEE=c(2,2),        # Máximo de erros esperados permitidos
                     truncQ=2,            # Corta se a qualidade cair abaixo de 2
                     rm.phix=TRUE,        # Remove genoma do vírus PhiX (controle comum)
                     compress=TRUE,       # Salva em formato compactado (.gz)
                     multithread=TRUE)    # Usa múltiplos núcleos do processador

# 5. Aprendizado das Taxas de Erro (Core do algoritmo DADA2)
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# 6. Desreplicação e Inferência de Amostras
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# 7. Salvar o progresso parcial para não perder o processamento pesado
save(dadaFs, dadaRs, errF, errR, file = "results/dada2_partial_output.RData")
print("Fase inicial do DADA2 concluída com sucesso!")
