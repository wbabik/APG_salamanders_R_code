# Here we build for each species ordinary linear multiple regression models
# with individual alpha diversities:
# APGs(and their various subsets) as response and MHC-I and average nonAPG
# as predictors
# we also plot APGs(MHC-I) with values for nonAPGs color-coded
# as input "*individual_alphadiv*" .rds files generated previously

library(tidyverse)

source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

#subsets sum_ind_div to get genes of interest for a model
#prepares df which has separate columns for PDs used as 
#response, explanatory and "covariate" variables
#this is ready for PGLS modeling
ind_div_to_xy <- function(arg_list){
  #gene as explanatory
  x <- get(arg_list[["x"]])
  #gene as response
  y <- get(arg_list[["y"]])
  #optional gene as "covariate"
  if (!is.null(arg_list[["z"]])) z <- get(arg_list[["z"]]) else z <- NULL
  #character vectors containing names of genes of interest
  xyz <- c(x, y, z)
  xz <- c(x, z)
  #get diversities for genes in xyz
  d <- sum_ind_div %>% ungroup() %>% filter(gene_sym %in% xyz)
  #get the response as separate df
  dy <- d %>% filter(gene_sym %in% y)
  #seed the result with response
  dxyz <- dy
  for (i in xz) {
    #add gene name for PDs variables
    alf <- sym(paste0("PD_alpha_", i))
    di <- d %>% filter(gene_sym %in% i) %>% select(sp_abr1, id, seq_type, PD_alpha) %>% 
      rename(!!alf := PD_alpha)
    #put everything together ready for modeling
    dxyz <- dxyz %>% left_join(di, by = c("sp_abr1", "id", "seq_type"))
  }
  return(dxyz)
}

fit_ols <- function (xyz, taxa) {
  d <- ind_div_to_xy(xyz)
  x <- get(xyz[["x"]])
  y <- get(xyz[["y"]])
  z <- if (!is.null(xyz[["z"]])) get(xyz[["z"]]) else NULL
  #print(z)
  #not optimal but the list is not very long
  res <- NULL
  for (t in taxa) {
    print(t)
    dt <- d %>% filter(sp_abr1 == t)
    #iterates through dependent variables (commonly only one) 
    for (yi in y) {
      print(yi)
      #iterates through sequence types, 
      #in this version DNA is left out
      for (st in c("codon", "AA", "AAGhm")) {
        print(st)
        #get the right subset of observations
        di <- dt %>% filter(gene_sym == yi, seq_type == st) %>% as.data.frame()
        if (nrow(di) < 15) {
          mi_res <- data.frame(species = t, form = NA, x = x, y = yi, cov = NA, 
                               seq_type = st, Py = NA, Pcov = NA,
                               y_est = NA, cov_est = NA, 
                               mod = NA)
        } else {
          dep <- "PD_alpha"
          e <- paste0("PD_alpha_", x)
          #two explanatory variables
          if (!is.null(z)) {
            e_cov <- paste0("PD_alpha_", z)
            f_main <- paste0(dep, " ~ ", e, " + ", e_cov)
            m_main <- lm(as.formula(f_main), data = di)
            t_table <- as.data.frame(summary(m_main)$coefficients)
            p <- t_table$`Pr(>|t|)`
            est <- t_table$`Estimate`
            mi_res <- data.frame(species = t, form = f_main, x = x, y = yi, cov = z, 
                                 seq_type = st, Py = p[[2]], Pcov = p[[3]],
                                 y_est = est[[2]], cov_est = est[[3]], 
                                 mod = I(list(main = m_main)))
            #one explanatory variable
          } else {
            f_main <- paste0(dep, " ~ ", e)
            m_main <- lm(as.formula(f_main), data = di)
            t_table <- as.data.frame(summary(m_main)$coefficients)
            p <- t_table$`Pr(>|t|)`
            est <- t_table$`Estimate`
            mi_res <- data.frame(species = t, form = f_main, x = x, y = yi, cov = NA, 
                                 seq_type = st, Py = p[[2]], Pcov = NA,
                                 y_est = est[[2]], cov_est = NA, 
                                 mod = I(list(main = m_main)))
          }
          
        }
        res <- rbind(res, mi_res)
      }
    }
  }
  return(res)
}

# BODY #

#Define constants ####
#names of the genes in subsets
APG <- c("PSMB8", "PSMB9", "TAP1", "TAP2", "TAPBP")
TAP <- c("TAP1", "TAP2")
PSMB <- c("PSMB8", "PSMB9")
APG_All <- c(APG, "APG", "PSMBs", "TAPs")
APG_subsets <- c("APG", "PSMBs", "TAPs", "TAPBP")
nonAPG <- c("BRD2", "DAXX", "KIFC1", "RGL2", "RXRBA")
APGmean <- c("APG")
nonAPGmean <- c("nonAPG")
MHC <- c("MHCI")
#this colors are currently not used, but may be handy at some time
g_colors <- c("PSMB8" = "darkorange", "PSMB9" = "darkorange4", "TAP1" = "green", "TAP2" = "green4", "TAPBP" = "mediumorchid1",
              "BRD2" = "gray60", "DAXX" = "turquoise1", "KIFC1" = "blue", "RGL2" = "gray30", "RXRBA" = "coral",
              "APG" = "blue", "nonAPG" = "red", "MHCI" = "magenta")

#Read data ####

#important: NaN values were present only for q > 0, not sure what's the reason,
#but we're using only q = 0, so it shouldn't be a problem
mip_ind <- readRDS("All_individual_diversities_excluding_stop_min20cov.rds")
mhc_brd_ind <- readRDS("All_sp_indivdual_alphadiv_MHD_BRD.rds")
taxonomy <- read_tsv("taxonomy.txt")
genes <- read.table("gene_ids.txt", sep = "\t", header = T, stringsAsFactors = F)
cds <- readRDS("CDS_lengths.rds")
brd_mhc_al <- readRDS("Funct_alleles_df_BRD2_MHC_ex2_ex3_all_genera.rds")

#Filter  and calculate per gene per individual diversities ####

#use only q = 0 and distances other than DNA
#remove "d" suffixes from some ids and code id as numeric
#join information on gene symbols and taxonomy
#MIPs
dind <- mip_ind %>% filter(seq_type != "DNA", q == 0) %>% select(-q) %>% 
  mutate(id = as.numeric(str_replace(id, "d", ""))) %>% 
  left_join(genes, by = c("gene")) %>% left_join(cds, by = c("species", "gene_sym"))
#MHC & BRD
mind <- mhc_brd_ind %>% filter (q == 0) %>% select(-q) %>% 
  rename(sp_abr1 = species) %>% left_join(taxonomy, by = "sp_abr1")

#per gene individual diversity
div_dind <- dind %>% group_by(species, id, class, gene_sym, seq_type) %>% 
  mutate(seg_weight = len/sum(len),
         alpha_weight = alpha*seg_weight) %>%
  summarise(reseq_len = sum(len),
            PD_alpha = sum(alpha_weight)) %>%
  rename(sp_abr1 = species) %>% left_join(taxonomy, by = "sp_abr1")

id_dind <- div_dind %>% pull(id) %>% unique()
#Per gene/group diversities ####

#arithmetic averages for PSMBs and TAPs
ind_psmb89 <- div_dind %>% ungroup() %>% filter(gene_sym %in% PSMB) %>% select(-gene_sym) %>%
  group_by(sp_abr1, class, seq_type, id, family, genus, species, gen_abr1, sp_abr2) %>%
  summarise(reseq_len = sum(reseq_len),
            PD_alpha = mean(PD_alpha, na.rm = TRUE)) %>%
  mutate(gene_sym = "PSMBs")

ind_tap12 <- div_dind %>% ungroup() %>% filter(gene_sym %in% TAP) %>% select(-gene_sym) %>%
  group_by(sp_abr1, class, seq_type, id, family, genus, species, gen_abr1, sp_abr2) %>%
  summarise(reseq_len = sum(reseq_len),
            PD_alpha = mean(PD_alpha, na.rm = TRUE)) %>%
  mutate(gene_sym = "TAPs")

ids_tap12 <- ind_tap12 %>% pull(id) %>% unique()
ids_psmb89 <- ind_psmb89 %>% pull(id) %>% unique()
#individuals that have PSMBs but not TAPs(some Triturus ind)
ids_to_exclude <- setdiff(ids_psmb89, ids_tap12)
#list of individuals that were successfully genotyped with MIPs
ids_mipped <- setdiff(id_dind, ids_to_exclude)


#combine the three above
#identical(sort(names(div_dind)), sort(names(ind_tap12)))
ind_divAPGnAPG <- bind_rows(div_dind, ind_psmb89, ind_tap12)

#process BRD2 and MHC ###
#for MHC involves weighted average across both exons
#get mean allele length for BRD2 and MHC exons
all_len <- brd_mhc_al %>% group_by(genus, species, exon) %>%
  summarise(mean_len_bp = mean(len, na.rm = TRUE))

ind_divBRD2MHC <- mind %>% filter(id %in% ids_mipped) %>% 
  left_join(all_len, by = c("genus", "species","exon")) %>% rename(PD_alpha = alpha)
to_del <- c("exon", "n_hap", "N_typed", "S", "S_unamb", "dmax")

ind_divBRD2 <- ind_divBRD2MHC %>% filter(exon == "brd2") %>% select(-all_of(to_del)) %>%
  add_column(class = "nonAPG", gene_sym = "BRD2", .after = "sp_abr1") %>%
  #relocate(mean_len_bp, PD_alpha, .after = "q") %>%
  rename(reseq_len = mean_len_bp)

##check whether names match
#identical(sort(names(ind_divAPGnAPG)), sort(names(ind_divBRD2)))

#consider including here only individuals with both MHC exons
#currently no such filter implemented
ind_divMHC <- ind_divBRD2MHC %>% filter(exon != "brd2") %>%
  select(-all_of(to_del))  %>%
  group_by(sp_abr1, family, genus, species, gen_abr1, sp_abr2, seq_type, id) %>%
  mutate(seg_weight = mean_len_bp/sum(mean_len_bp),
         PD_alpha_weight = PD_alpha*seg_weight) %>%
  summarise(reseq_len = sum(mean_len_bp),
            PD_alpha = sum(PD_alpha_weight)) %>%
  mutate(gene_sym = "MHCI",
         class = "MHC")

sum_ind_div <- bind_rows(ind_divAPGnAPG, ind_divMHC, ind_divBRD2) %>% select(-c(gen_abr1, sp_abr2)) %>% ungroup()

#calculates average diversities for APGs and nonAPGs
temp <- sum_ind_div %>% group_by(sp_abr1, class, seq_type, family, genus, species, id) %>%
  summarise(reseq_len = sum(reseq_len),
            PD_alpha = mean(PD_alpha, na.rm = TRUE)) %>%
  filter(class != "MHC") %>%
  mutate(gene_sym = class)

#adds average diversities calculated above to the diversity table
#this is the complete diversity table
#can be extended if needed, e.g. with weighted averages
sum_ind_div <- bind_rows(sum_ind_div, temp) %>% 
  filter(id %in% ids_mipped) %>% mutate(taxon = paste0(genus, "_", species))

#Model ####

#list of taxa to include
taxa <- sum_ind_div %>% pull(sp_abr1) %>% unique()

#list of gene sets to be used in modeling
geny <- list(x = "MHC", y = "APG_subsets", z = "nonAPGmean")


#ble <- ind_div_to_xy(geny)
#model, the resulting df contains also the model itself
#NA are things where modeling was not possible because of no data (for example gene loss)
species_alpha_models <- fit_ols(geny, taxa) %>% filter(!is.na(form))
#as above but model removed from df
sam <- species_alpha_models %>% select(-mod)

#write_lnx_head(select(species_alpha_models, -mod), "Alpha_diversities_within_species_modeling.txt")


#Prepare for plotting ####
#loops through model results df
#checks if MHC-I was significant
#if so exctracts regression line parameters
#extracts raw data from model
#joins with metainfo
#the result is dataframe ready to plot
ind_to_fig_plot <- vector("list", nrow(species_alpha_models))
for (i in 1:nrow(species_alpha_models)) {
  r <- species_alpha_models[i, ]
  m <- r$mod[[1]]
  print(summary(m))
  pval <- summary(m)$coefficients[2,4]
  print(pval)
  #ax +b
  #if a not significant, then NA
  if (is.nan(pval) || pval > 0.05) {
    a <- NA
    b <- NA
  } else {
    a <- m$coefficients[[2]]
    b <- m$coefficients[[1]]
  }
  print(a)
  print(b)
  if (!is.na(m)){
    d <- as.data.frame(m$model)
    #add NA if not significant
    data <- data.frame(response = d$PD_alpha, MHC_I = d$PD_alpha_MHCI, nonAPG = d$PD_alpha_nonAPG, intercept = b, slope = a)
    i_df <- cbind(select(r, -mod), data)
    ind_to_fig_plot[[i]] <- i_df
  }
}

ind_to_fig_plot <- bind_rows(ind_to_fig_plot)


#Plot ####
theme_set(theme_bw(base_size = 14) + theme(plot.title = element_text(hjust = 0.5, size = 20),
                                           axis.title = element_text(size = 16),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank())) 
update_geom_defaults("point", list(size = 1))
update_geom_defaults("abline", list(size = 1, colour = "grey90"))

AA <- ind_to_fig_plot %>% filter(seq_type == "AA")


fAA <- ggplot(AA, aes(x = MHC_I, y = response, colour = y)) + geom_point() + 
  geom_abline(aes(intercept = intercept, slope=slope, colour = y)) + 
  facet_wrap(vars(species), scales = "free") +
  labs(title = "Individual-level diversity APG vs MHC-I (protein p-distance)") +
  xlab("MHC-I") + ylab("APG") + theme(axis.title.x = element_text(face = "italic"))

fAA
ggsave("Individual_APG_vs_MHC_diversity.pdf", width = 15, height = 12, units = "in")
