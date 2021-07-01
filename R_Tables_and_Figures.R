#Load libraries and functions####

library(tidyverse)

#load common functions
source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

#Load data ####
ind_inf <- read_tsv("individual_info.txt")
localities <- read_csv("localities.csv")
seg_div <- readRDS("All_sp_div_excluding_stop_min20cov.rds")
div_sum <- readRDS("Diversities_summary_long_format.rds")

cds <- readRDS("CDS_lengths.rds")
genes <- read_tsv("gene_ids.txt") 
taxonomy <- read_tsv("taxonomy.txt") %>% select(-c(gen_abr1, sp_abr2))

#Define outputs ####
sum_sp_g_file <- "Resequencing_summary_APG_nonAPG_per_gene_per_species.txt"
sum_g_file <- "Resequencing_summary_APG_nonAPG_per_gene.txt"

#Data to summarise ####
d <- seg_div %>% left_join(genes, by = c("gene")) %>% left_join(cds, by = c("species", "gene_sym")) %>% 
  filter(q == 0 & (seq_type == "DNA" | seq_type == "codon")) %>% 
  select(species, gene_sym, class, cds_length, seq_type, segment, len, cds_len_bp, n_hap, N_typed, fr_typed, S, S_unamb) %>% 
  mutate(cds_len_bp = as.numeric(cds_len_bp))

#Summaries per gene per species ####
sum_sp_g <- d %>% filter(fr_typed >= 0.5) %>%  group_by(species, class, gene_sym, cds_length) %>% 
  summarise(cds_length = unique(cds_length),
            mean_seg_len = mean(len[seq_type == "codon"], na.rm = TRUE),
            median_seg_len = median(len[seq_type == "DNA"],na.rm = TRUE),
            mean_seg_codon_bp_len = 3*mean(len[seq_type == "codon"], na.rm = TRUE),
            median_seg_codon_bp_len = 3*median(len[seq_type == "codon"], na.rm = TRUE),
            cds_reseq = sum(cds_len_bp[seq_type == "DNA"]),
            codon_bp_reseq = 3*sum(len[seq_type == "codon"]),
            non_cds_reseq = sum(len[seq_type == "DNA"]) - cds_reseq,
            fr_cds_reseq = cds_reseq/unique(cds_length),
            fr_codon_reseq = codon_bp_reseq/unique(cds_length))

c_to_sum_g <- c("cds_reseq", "codon_bp_reseq", "fr_cds_reseq", "fr_codon_reseq")

fun_to_sum <- list(min = ~min(.x),
                   mean = ~mean(.x),
                   sd = ~sd(.x),
                   median = ~median(.x),
                   max = ~max(.x))

#Summaries per species ####
sum_g <- sum_sp_g %>% group_by(class, gene_sym) %>% 
  summarise(n_spec = n(),
            mean_length = mean(cds_length),
            across(all_of(c_to_sum_g), fun_to_sum))

#Write summaries to text ####
write_lnx_head(sum_sp_g, sum_sp_g_file)
write_lnx_head(sum_g, sum_g_file)


#Table 1: per species per class diversities 
#for alpha q = 0, for gamma q = 1 ####

div_sum_temp <- div_sum %>%
  select(family, genus, species, seq_type, gene_sym, q, PD_alpha, PD_gamma) %>%
  filter(q %in% c(0, 1), gene_sym %in% c("APG", "nonAPG", "MHCI"), seq_type != "DNA") 
#%>% select(-q)

table_1 <- div_sum_temp %>% pivot_wider(names_from = c(seq_type, gene_sym, q), values_from = c(PD_alpha, PD_gamma)) %>%
  select(family, genus, species, order(colnames(.))) %>% 
  select(-((starts_with("PD_alpha") & ends_with("1")) | (starts_with("PD_gamma") & ends_with("0")))) 
write_tsv(table_1, "Table_1.txt")


#Generate several supplementary tables ####

#Table S1 - samples 
ind_loc <- ind_inf %>% select(-country) %>% left_join(localities, by = "locality") %>% 
  left_join(taxonomy, by = c("genus", "species"))
mipped <- readRDS("All_individual_diversities_excluding_stop_min20cov.rds") %>% pull(id) %>% unique() %>% 
  str_replace("d", "") %>% as.numeric() %>% as.data.frame()
names(mipped) <- "id"
table_s1 <- mipped %>% left_join(ind_loc, by = "id") %>% select(-subspecies) %>% 
  group_by(family, genus, species, locality, country, latitude, n_s, longitude, e_w) %>% summarise(N = n())
write_lnx_head(table_s1, "Table_S1.txt")

#BRD2 assumed mean cds length 2500 bp (there are slight differences between taxa)
#resequenced length in 22 species is 300 bp and in 8 297
BRD2dat <- c(rep(300, 22), rep(297, 8))
meanBRD <- mean(BRD2dat)
sdBRD <- sd(BRD2dat)
lenBRD <- 2500.0

spBRDdf <- data.frame(species = unique(sum_sp_g$species), class = "nonAPG", gene_sym = "BRD2", cds_length = 2500) %>% 
  mutate(cds_reseq = ifelse(grepl("Tri_|Sal", species), 297, 300))

temp <- sum_sp_g %>% select(species, class, gene_sym, cds_length, cds_reseq) %>% 
  bind_rows(spBRDdf) %>% rename(sp_abr1 = species) %>% left_join(taxonomy, by = "sp_abr1")

#Table S4 resequencing stats
table_S4 <- temp %>% group_by(family, genus, species, class) %>% 
  summarise(`resequenced [bp]` = sum(cds_reseq),
  `fraction cds resequenced` = sum(cds_reseq)/sum(cds_length)) %>% 
  pivot_wider(names_from = class, values_from = c(`resequenced [bp]`, `fraction cds resequenced`, )) %>% 
  select(family:species, 4, 6, 5, 7) %>% 
  rename()
write_tsv(table_S4, "Table_S4.txt")

#Table S5 per gene summary
table_S5 <- sum_g %>% ungroup() %>% select(class, gene_sym, n_spec, mean_length, cds_reseq_mean, cds_reseq_sd) %>% 
  add_row(class = "nonAPG", gene_sym = "BRD2", n_spec = 30.0, 
          mean_length = lenBRD, cds_reseq_mean = meanBRD, cds_reseq_sd = sdBRD) %>% 
  arrange(class, gene_sym) %>% 
  rename(category = class, gene = gene_sym, `N species` = n_spec, 
         `Mean cds length [bp]` = mean_length, `Mean cds resequenced [bp]` = cds_reseq_mean, 
         `SD cds resequenced [bp]` = cds_reseq_sd, )
write_tsv(table_S5, "Table_S5.txt")


#Table S6 per gene stats
table_S6 <- temp %>% ungroup() %>% select(-sp_abr1) %>% 
  mutate(`fraction resequenced` = cds_reseq/cds_length) %>% 
  select(family:species, class:`fraction resequenced`) %>% 
  rename(category = class, gene = gene_sym, `cds length [bp]` = cds_length, `resequenced length [bp]` = cds_reseq) %>% 
  arrange(category, gene, family, genus, species)

write_tsv(table_S6, "Table_S6.txt")

#Table S8 per gene diversity
table_S8 <- div_sum %>% select(-c(sp_abr1, reseq_len, taxon)) %>% select(family:species, class:PD_gamma) %>%
  rename(category = class, gene = gene_sym, distance = seq_type, 
         `within-individual (alpha) dversity`= PD_alpha,
         `species-wide (gamma) diversity` = PD_gamma) %>% 
  filter(distance != "DNA", category !="ex3") %>% mutate(distance = ifelse(distance=="codon", "dS", distance)) %>% 
  arrange(category, gene, distance, q, family, genus, species)
  
write_tsv(table_S8, "Table_S8.txt")


#Tables S9-S11 modeling
s9 <- read_tsv("model_characteristics_PGLS_modeling_results.txt") %>% mutate(table = "s9")
s10 <- read_tsv("model_characteristics_PGLS_modeling_results_15_best_covered.txt") %>% mutate(table = "s10")
s11 <- read_tsv("model_characteristics_PGLS_modeling_results_classical.txt") %>%  mutate(table = "s11")
ss <- bind_rows(s9, s10, s11) %>%   filter(diversity == "gamma" | (diversity == "alpha" & q == 0)) %>% 
  filter(predictor != "ex3", response != "ex3") %>%
  arrange(table, predictor, response, diversity, distance, q, parameter)

ss %>% filter(table == "s9") %>% select(-table) %>% write_tsv("Table_S9.txt")

ss %>% filter(table == "s10") %>% select(-table) %>% write_tsv("Table_S10.txt")

ss %>% filter(table == "s11") %>% select(-table) %>% write_tsv("Table_S11.txt")
