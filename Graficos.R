###########################################################################################
#                                Pacotes necessários                                      #
###########################################################################################

library(phyloseq)     # Manipulação de dados de microbioma
library(ggplot2)      # Gráficos estáticos
library(dplyr)        # Manipulação de dados
library(tidyr)        # Transformações tidyr (separar, unir colunas)
library(readr)        # Leitura de ficheiros .tsv
library(tibble)       # Manipulação de data frames como tibbles
library(vegan)        # Análise de bray-curtis
library(plotly)       # Gráficos interativos em 3D
library(ggpubr)       # Boxplot + testes estatísticos

###########################################################################################
#                                Leitura e preparação dos dados                           #
###########################################################################################

# Leitura dos ficheiros de taxonomia e contagens
taxonomy <- read_tsv("/caminho/para/OTUs.97.cons.taxonomy")
taxonomy_table <- read_tsv("/caminho/para/OTUs.97.rep.count_table")

# Limpeza da taxonomia e separação dos níveis taxonómicos
taxonomy$Taxonomy <- gsub(";$", "", taxonomy$Taxonomy)
tax_split <- taxonomy %>%
  separate(Taxonomy, into = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"),
           sep = ";", fill = "right", extra = "drop") %>%
  column_to_rownames("OTU")
tax_matrix <- as.matrix(tax_split)

# Processar tabela de contagens
taxonomy_table <- taxonomy_table %>%
  column_to_rownames("Representative_Sequence")
counts_matrix <- as.matrix(taxonomy_table)
mode(counts_matrix) <- "numeric"

# Criar objeto phyloseq inicial
ps <- phyloseq(
  otu_table(counts_matrix, taxa_are_rows = TRUE),
  tax_table(tax_matrix)
)

###########################################################################################
#                            Leitura e preparação dos metadados                           #
###########################################################################################

sra_metadata <- read.csv("/caminho/para/UTF-8assigment2.csv", stringsAsFactors = FALSE)
rownames(sra_metadata) <- sra_metadata$host_subject_id
counts_matrix <- counts_matrix[, -1]  # Remove coluna 'total'

# Interseção das amostras entre contagens e metadados
common_samples <- intersect(colnames(counts_matrix), rownames(sra_metadata))
sra_metadata_filtered <- sra_metadata[match(common_samples, rownames(sra_metadata)), , drop = FALSE]
stopifnot(nrow(sra_metadata_filtered) > 0)

# Adicionar metadados ao objeto phyloseq
sample_data(ps) <- sample_data(sra_metadata_filtered)
sample_data(ps)$Group <- factor(sample_data(ps)$Host_disease)

# Filtrar grupos de interesse
ps_filtered <- subset_samples(ps, Group %in% c("Negative", "Blank"))

###########################################################################################
#                      Diversidade Alfa + Testes Estatísticos                             #
###########################################################################################

# Preparar metadados
meta_df <- as(sample_data(ps_filtered), "data.frame") %>%
  rownames_to_column("SampleID")

# Leitura dos dados de diversidade alfa
aalpha_div <- read_tsv("/caminho/para/alpha-div.tsv") %>%
  rename(SampleID = group)

# Junção dos dados de diversidade com os metadados
alpha_div_joined <- inner_join(aalpha_div, meta_df, by = "SampleID")

# Teste de normalidade Shapiro-Wilk
shannon_blank <- alpha_div_joined$shannon[alpha_div_joined$Group == "Blank"]
shannon_neg <- alpha_div_joined$shannon[alpha_div_joined$Group == "Negative"]
shapiro_blank <- if(length(shannon_blank) >= 3) shapiro.test(shannon_blank) else NULL
shapiro_neg <- if(length(shannon_neg) >= 3) shapiro.test(shannon_neg) else NULL

cat("Teste Shapiro-Wilk para Blank:\n")
if (!is.null(shapiro_blank)) print(shapiro_blank) else cat("Insuficiente.\n")
cat("Teste Shapiro-Wilk para Negative:\n")
if (!is.null(shapiro_neg)) print(shapiro_neg) else cat("Insuficiente.\n")

# Teste de comparação entre grupos
if (!is.null(shapiro_blank) && !is.null(shapiro_neg) &&
    shapiro_blank$p.value > 0.05 && shapiro_neg$p.value > 0.05) {
  test_result <- t.test(shannon ~ Group, data = alpha_div_joined)
} else {
  test_result <- wilcox.test(shannon ~ Group, data = alpha_div_joined)
}
print(test_result)

###########################################################################################
#                                     Gráfico A (Beta-diversidade)                        #
###########################################################################################

beta_div <- read_tsv("/caminho/para/beta-div.tsv") 

# Extrair eigenvalues
ordination$values$Eigenvalues

# Calcular a ordinação com método PCoA e distância de Bray-Curtis
ordination <- ordinate(ps_filtered, method = "PCoA", distance = "bray")

# Extrair coordenadas (Axis.1, Axis.2, Axis.3)
ordination_df <- as.data.frame(ordination$vectors) %>%
  rownames_to_column("SampleID") %>%
  left_join(meta_df, by = "SampleID")
colnames(ordination_df)
head(ordination_df)

# Definir cores e formas para os grupos
group_colors <- c("Negative" = "#000", "Blank" = "#FF007F")
group_shapes <- c("Negative" = "circle", "Blank" = "circle")

# Gráfico interativo 3D com Plotly
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
    zaxis = list(title = "Axis 3")))

###########################################################################################
#                                     Gráfico B (Phylum)                                  #
###########################################################################################

ps_grouped_phylum <- tax_glom(ps_filtered, taxrank = "Phylum")
phylum_names <- as.character(tax_table(ps_grouped_phylum)[, "Phylum"])
phylum_names[is.na(phylum_names) | phylum_names == ""] <- "Unclassified"
tax_table(ps_grouped_phylum)[, "Phylum"] <- phylum_names

ps_phylum_relabund <- transform_sample_counts(ps_grouped_phylum, function(x) x / sum(x))
df_phylum <- psmelt(ps_phylum_relabund)

df_phylum_mean <- df_phylum %>%
  group_by(Group, Phylum) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

ggplot(df_phylum_mean, aes(x = Group, y = Abundance * 100, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Grupo", y = "Abundância Relativa Média (%)", fill = "Phylum") +
  theme_minimal()

###########################################################################################
#                                     Gráfico C (Genus)                                   #
###########################################################################################

ps_grouped_rel <- transform_sample_counts(ps_grouped, function(x) x / sum(x) * 100)
ps_grouped_genus <- tax_glom(ps_grouped_rel, taxrank = "Genus")

genus_names <- as.character(tax_table(ps_grouped_genus)[, "Genus"])
genus_names[is.na(genus_names) | genus_names == ""] <- "Unclassified"
tax_table(ps_grouped_genus)[, "Genus"] <- genus_names

df <- psmelt(ps_grouped_genus)
df_mean <- df %>%
  group_by(Group, Genus) %>%
  summarise(Abundance = mean(Abundance), .groups = "drop")

cutoff <- 0.01
genera_keep <- df_mean %>%
  group_by(Genus) %>%
  summarise(Total = sum(Abundance)) %>%
  filter(Total >= cutoff) %>%
  pull(Genus)

df_final <- df_mean %>%
  mutate(Genus = ifelse(Genus %in% genera_keep, Genus, "Other")) %>%
  group_by(Group, Genus) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

ggplot(df_final, aes(x = Group, y = Abundance * 100, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Grupo", y = "Abundância Relativa Média (%)", fill = "Genus") +
  theme_minimal()

###########################################################################################
#                                     Gráfico D (Razão P/A)                               #
###########################################################################################

ps_phylum <- tax_glom(ps, taxrank = "Phylum")
ps_rel <- transform_sample_counts(ps_phylum, function(x) x / sum(x))
df <- psmelt(ps_rel)

df_sub <- df %>%
  filter(Phylum %in% c("Proteobacteria", "Actinobacteriota")) %>%
  select(Sample, Phylum, Abundance, Group)

df_summarised <- df_sub %>%
  group_by(Sample, Phylum, Group) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop")

df_wide <- df_summarised %>%
  pivot_wider(names_from = Phylum, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  mutate(Ratio = Proteobacteria / (Actinobacteriota + 1e-6))

ggplot(df_wide, aes(x = Group, y = Ratio, fill = Group)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.6) +
  labs(x = "Grupo", y = "Razão Proteobacteria / Actinobacteria") +
  theme_minimal()

###########################################################################################
#                                     Gráfico E (Géneros ANCOM)                           #
###########################################################################################

ancom_genera <- c("Streptococcus", "Corynebacterium", "Anaerococcus", 
                  "Finegoldia", "Enterococcus", "Bacillus")

ps_rel <- transform_sample_counts(ps_grouped_genus, function(x) x / sum(x) * 100)
genus_names <- as.character(tax_table(ps_rel)[, "Genus"])
genus_names[is.na(genus_names) | genus_names == ""] <- "Unclassified"
tax_table(ps_rel)[, "Genus"] <- genus_names

keep_taxa <- taxa_names(ps_rel)[tax_table(ps_rel)[, "Genus"] %in% ancom_genera]
ps_ancom <- prune_taxa(keep_taxa, ps_rel)

df <- psmelt(ps_ancom) %>%
  group_by(Genus, Group) %>%
  summarise(
    Mean_Abundance = mean(Abundance),
    SD = sd(Abundance),
    .groups = "drop"
  )

ggplot(df, aes(x = Group, y = Mean_Abundance, fill = Group)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.6, alpha = 0.8) +
  geom_errorbar(aes(ymin = Mean_Abundance - SD, ymax = Mean_Abundance + SD),
                width = 0.2, position = position_dodge(0.6)) +
  facet_wrap(~ Genus, scales = "free_y", ncol = 3) +
  labs(x = "", y = "Mean Relative Abundance (%)") +
  scale_fill_manual(values = c("Negative" = "#000000", "Blank" = "#FF007F")) +
  theme_minimal()
