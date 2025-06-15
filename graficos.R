# Pacotes necessários
library(phyloseq)     # Manipulação de dados de microbioma
library(ggplot2)      # Gráficos estáticos
library(dplyr)        # Manipulação de dados
library(tidyr)        # Transformações tidyr (separar, unir colunas)
library(readr)        # Leitura de ficheiros .tsv
library(tibble)       # Para manipular data frames como tibbles
library(vegan)        # Análise de ecologia (ex: distâncias Bray-Curtis)
library(plotly)       # Gráficos interativos em 3D
library(ggpubr)       # Boxplot + testes estatísticos

# Carregar dados de taxonomia e contagens
taxonomy <- read_tsv("/home/sara-sampaio/16S-pipeline_outputs/results/consensus/OTUs.97.cons.taxonomy")  
taxonomy_table <- read_tsv("/home/sara-sampaio/16S-pipeline_outputs/results/consensus/OTUs.97.rep.count_table")  

# Remover ponto e vírgula no final da string de taxonomia
taxonomy$Taxonomy <- gsub(";$", "", taxonomy$Taxonomy)

# Separar taxonomia por níveis e definir OTU como rowname
tax_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";", fill = "right", extra = "drop") %>%
  column_to_rownames("OTU")

# Converter tabela de taxonomia para matriz (necessária para phyloseq)
tax_matrix <- as.matrix(tax_split)

# Processar tabela de contagens e converter para matriz
taxonomy_table <- taxonomy_table %>%
  column_to_rownames("Representative_Sequence")
counts_matrix <- as.matrix(taxonomy_table)
mode(counts_matrix) <- "numeric"  # Garantir que são valores numéricos

# Criar objeto phyloseq com matriz de contagens e taxonomia
ps <- phyloseq(
  otu_table(counts_matrix, taxa_are_rows = TRUE),
  tax_table(tax_matrix)
)

# Criar metadados simples com base no nome da amostra
sample_data_df <- data.frame(
  Group = factor(
    ifelse(grepl("Neg", colnames(counts_matrix), ignore.case = TRUE), "Negative",
           ifelse(grepl("Blank", colnames(counts_matrix), ignore.case = TRUE), "Blank", "Other")),
    levels = c("Negative", "Blank", "Other")  # <- isto garante que todos os níveis existem
  ),
  row.names = colnames(counts_matrix)
)
sample_data_df$Group <- factor(sample_data_df$Group)

# Só agora adicionar os metadados ao objeto phyloseq
sample_data(ps) <- sample_data(sample_data_df)

# Filtrar apenas grupos de interesse
ps_filtered <- subset_samples(ps, Group %in% c("Negative", "Blank"))
table(sample_data(ps_filtered))
sample_names(ps_filtered)

###########################################################################################
#                      Diversidade Alfa + Testes Estatísticos                             #
###########################################################################################

#  Preparar metadados para amostras filtradas (ps_filtered) 
meta_df <- as(sample_data(ps_filtered), "data.frame") %>%   # converter para data.frame
  tibble::rownames_to_column("SampleID")                    # transformar rownames em coluna

# Verificar as amostras e grupos
print(meta_df$SampleID)
print(table(meta_df$Group))

# Carregar o ficheiro alfa
alpha_div <- readr::read_tsv("/home/sara-sampaio/16S-pipeline_outputs/postprocessing/alpha-div.tsv") %>%
  dplyr::rename(SampleID = group)
sample_names(ps_filtered)
# Verificar amostras em alpha_div
print(unique(alpha_div$SampleID))

# --- Fazer inner join entre alpha-diversity e metadados ---
alpha_div_joined <- dplyr::inner_join(alpha_div, meta_df, by = "SampleID")

# Confirmar número de amostras por grupo depois do join
print(table(alpha_div_joined$Group))

# --- Testes estatísticos sobre índice Shannon ---

# Separar dados por grupo para Shapiro-Wilk
shannon_blank <- alpha_div_joined$shannon[alpha_div_joined$Group == "Blank"]
shannon_neg <- alpha_div_joined$shannon[alpha_div_joined$Group == "Negative"]

# Teste de normalidade só se tiver 3 ou mais amostras
shapiro_blank <- if(length(shannon_blank) >= 3) shapiro.test(shannon_blank) else NULL
shapiro_neg <- if(length(shannon_neg) >= 3) shapiro.test(shannon_neg) else NULL

cat("Teste Shapiro-Wilk para Blank:\n")
if(!is.null(shapiro_blank)) {
  print(shapiro_blank)
} else {
  cat("Número insuficiente de amostras para teste Shapiro-Wilk (Blank).\n")
}

cat("\nTeste Shapiro-Wilk para Negative:\n")
if(!is.null(shapiro_neg)) {
  print(shapiro_neg)
} else {
  cat("Número insuficiente de amostras para teste Shapiro-Wilk (Negative).\n")
}

# Escolher teste estatístico conforme normalidade e disponibilidade de teste
if(!is.null(shapiro_blank) && !is.null(shapiro_neg) && 
   shapiro_blank$p.value > 0.05 && shapiro_neg$p.value > 0.05) {
  
  cat("\nAmbos os grupos seguem distribuição normal. Usando teste t não pareado.\n")
  test_result <- t.test(shannon ~ Group, data = alpha_div_joined)
  
} else {
  
  cat("\nPelo menos um grupo não tem normalidade ou não há amostras suficientes. Usando teste de Wilcoxon.\n")
  test_result <- wilcox.test(shannon ~ Group, data = alpha_div_joined)
  
}

cat("\nResultado do teste de comparação entre grupos:\n")
print(test_result)

###########################################################################################
#                                     Gráfico A                                           #
###########################################################################################

# Calcular matriz de distância
distance_matrix <- phyloseq::distance(ps_filtered, method = "bray")  
ordination <- ordinate(ps, method = "PCoA", distance = distance_matrix) 

# Extrair coordenadas
ordination_df <- as.data.frame(ordination$vectors)  
ordination_df$SampleID <- rownames(ordination_df)

# Juntar coordenadas com metadados
meta <- as(sample_data(ps_filtered), "data.frame")
meta$SampleID <- rownames(meta)
ordination_df <- left_join(ordination_df, meta, by = "SampleID")

# Definir cores e formas para o gráfico
group_colors <- c("Negative" = "#000", "Blank" = "#FF007F")  
group_shapes <- c("Negative" = "circle", "Blank" = "circle")  

# Criar gráfico interativo em 3D
plot_ly(ordination_df, 
        x = ~Axis.1, y = ~Axis.2, z = ~Axis.3,
        color = ~Group, colors = group_colors,
        symbol = ~Group, symbols = group_shapes,
        text = ~paste("Sample:", SampleID, "<br>Group:", Group),
        hoverinfo = 'text') %>%
  add_markers(size = 10) %>%
  layout(scene = list(
    xaxis = list(title = "Axis 1"),
    yaxis = list(title = "Axis 2"),
    zaxis = list(title = "Axis 3"))
  )

###########################################################################################
#                                         Gráfico B                                       #
###########################################################################################

# Agrupar por phylum
ps_grouped_phylum <- tax_glom(ps_filtered, taxrank = "Phylum")  

# Substituir valores NA por "Unclassified"
phylum_names <- as.character(tax_table(ps_grouped_phylum)[, "Phylum"])
phylum_names[is.na(phylum_names) | phylum_names == ""] <- "Unclassified"
tax_table(ps_grouped_phylum)[, "Phylum"] <- phylum_names

# Converter para abundância relativa
ps_phylum_relabund <- transform_sample_counts(ps_grouped_phylum, function(x) x / sum(x))
df_phylum <- psmelt(ps_phylum_relabund)

# Calcular média por grupo e phylum
df_phylum_mean <- df_phylum %>%
  group_by(Group, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

# Gráfico de barras empilhadas
ggplot(df_phylum_mean, aes(x = Group, y = Abundance * 100, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Grupo", y = "Abundância Relativa Média (%)", fill = "Phylum") +
  theme_minimal()

###########################################################################################
#                                        Gráfico C                                        #
###########################################################################################
ps_grouped_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x) * 100)
ps_grouped_genus <- tax_glom(ps_grouped_rel, taxrank = "Genus")

# Corrigir nomes em falta
genus_names <- as.character(tax_table(ps_grouped_genus)[, "Genus"])
genus_names[is.na(genus_names) | genus_names == ""] <- "Unclassified"
tax_table(ps_grouped_genus)[, "Genus"] <- genus_names

# Calcular abundância relativa média por género
ps_relabund <- transform_sample_counts(ps_grouped_genus, function(x) x / sum(x))
df <- psmelt(ps_relabund)

df_mean <- df %>%
  group_by(Group, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

# Filtrar géneros com abundância total < 1%
cutoff <- 0.01
genera_keep <- df_mean %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  filter(Total >= cutoff) %>%
  pull(Genus)

# Colapsar géneros "menos abundantes" como "Other"
df_final <- df_mean %>%
  mutate(Genus = ifelse(Genus %in% genera_keep, Genus, "Other")) %>%
  group_by(Group, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

ggplot(df_final, aes(x = Group, y = Abundance * 100, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Grupo", y = "Abundância Relativa Média (%)", fill = "Genus") +
  theme_minimal()

###########################################################################################
#                                         Gráfico D                                       #
###########################################################################################
ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
df <- psmelt(ps_rel)

# Selecionar apenas os dois phyla de interesse
df_sub <- df %>%
  filter(Phylum %in% c("Proteobacteria", "Actinobacteriota")) %>%
  select(Sample, Phylum, Abundance, Group)

# Calcular abundância por amostra
df_summarised <- df_sub %>%
  group_by(Sample, Phylum, Group) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

# Pivot para formato largo e calcular razão
df_wide <- df_summarised %>%
  pivot_wider(names_from = Phylum, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  mutate(Ratio = Proteobacteria / (Actinobacteriota + 1e-6))

# Boxplot 
ggplot(df_wide, aes(x = Group, y = Ratio, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(x = "Grupo", y = "Razão Proteobacteria / Actinobacteria") +
  theme_minimal()

###########################################################################################
#                                         Gráfico E                                       #
###########################################################################################
ancom_genera <- c("Streptococcus", "Corynebacterium", "Anaerococcus", 
                  "Finegoldia", "Enterococcus", "Bacillus")

# Converter para abundância relativa (%)
ps_rel <- transform_sample_counts(ps_grouped_genus, function(x) x / sum(x) * 100)
genus_names <- as.character(tax_table(ps_rel)[, "Genus"])
genus_names[is.na(genus_names) | genus_names == ""] <- "Unclassified"
tax_table(ps_rel)[, "Genus"] <- genus_names

# Selecionar apenas os géneros de interesse
keep_taxa <- taxa_names(ps_rel)[tax_table(ps_rel)[, "Genus"] %in% ancom_genera]
ps_ancom <- prune_taxa(keep_taxa, ps_rel)

# Calcular média e desvio padrão por grupo e género
df <- psmelt(ps_ancom) %>%
  group_by(Genus, Group) %>%
  summarise(
    Mean_Abundance = mean(Abundance),
    SD = sd(Abundance),
    .groups = "drop"
  )

# Gráfico com barras e erro padrão
ggplot(df, aes(x = Group, y = Mean_Abundance, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Abundance - SD, ymax = Mean_Abundance + SD),
                width = 0.2, position = position_dodge(0.6)) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
  labs(x = "", y = "Mean Relative Abundance (%)") +
  scale_fill_manual(values = c("Negative" = "#000000", "Blank" = "#FF007F")) +
  theme_minimal()
