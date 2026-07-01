# Pipeline de Análise de Diversidade Microbiana via Metabarcoding (16S/ITS) 🧬💻

Este repositório contém um fluxo de trabalho (*pipeline*) estruturado e automatizado para o processamento de dados brutos de sequenciamento de amplicons (Metabarcoding), cobrindo desde o controlo de qualidade das leituras até a modelação ecológica e visualização de comunidades microbianas.

O foco deste projeto é demonstrar a aplicação prática de ferramentas de bioinformática e ciência de dados no estudo de ecologia microbiana, consórcios biológicos e genómica ambiental.

---

## 📂 Estrutura do Repositório

O projeto segue as boas práticas de organização de diretórios em bioinformática:

* `data/`
    * `raw_data/`: Diretório destinado aos arquivos brutos de sequenciamento (`.fastq` ou `.fastq.gz`).
    * `reference_db/`: Local para armazenamento dos bancos de dados taxonómicos de referência (ex: SILVA, UNITE).
* `scripts/`: Scripts em linguagem R para execução modular do pipeline.
* `results/`
    * `figures/`: Gráficos exportados de diversidade alfa, beta e composição taxonómica.

---

## ⚙️ Fluxo de Trabalho e Módulos

O pipeline está dividido em três etapas principais executadas de forma sequencial:

### 1. Filtragem e Inferência de Variantes (DADA2)
O script `01_dada2_processing.R` implementa o algoritmo do pacote **DADA2**. Diferente dos métodos tradicionais baseados em OTUs (unidades taxonómicas operacionais), este pipeline trabalha com **ASVs (Amplicon Sequence Variants)**, garantindo a resolução de um único nucleotídeo através de:
* Controlo de qualidade e corte de reads por score de qualidade (`filterAndTrim`).
* Modelagem estatística adaptativa das taxas de erro da corrida de sequenciamento (`learnErrors`).
* Desreplicação e inferência exata de sequências biológicas reais, discriminando erros de leitura do sequenciador.

### 2. Resolução Taxonómica
O script `02_dada2_taxonomy.R` realiza a consolidação das matrizes biológicas:
* Alinhamento e junção das leituras pareadas (*Forward* e *Reverse*).
* Remoção *de novo* de quimeras (artefatos de PCR).
* Atribuição taxonómica utilizando algoritmos de classificação baseados em bancos de dados de referência (ex: taxonomia de bactérias ou fungos).

### 3. Ecologia Microbiana e Análise Estatística Downstream
O script `03_microbial_ecology.R` integra as matrizes geradas utilizando os pacotes **Phyloseq**, **Vegan** e **GGPlot2** para extrair métricas ecológicas:
* **Diversidade Alfa:** Cálculo de índices de riqueza observada e Shannon para avaliar a complexidade interna das amostras.
* **Diversidade Beta:** Análises de ordenação espacial (NMDS utilizando distância de Bray-Curtis) para avaliar a variação na composição das comunidades entre diferentes tratamentos ou localizações.
* **Perfil Taxonómico:** Geração de gráficos de abundância relativa acumulada para visualização da estrutura da comunidade a nível de Filo/Família.

---

## 🛠️ Tecnologias e Ferramentas Utilizadas

* **Linguagem Principal:** R (v4.x)
* **Ecossistema Ómico:** BioConductor (`dada2`, `phyloseq`)
* **Análise Estatística e Gráficos:** `vegan`, `ggplot2`
* **Controlo de Versão:** Git e GitHub

---
