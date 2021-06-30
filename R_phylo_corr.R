library(ape)
library(caper)
library(ggpubr)
library(tidyverse)

#load common functions
source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

# CONSTANTS ####

#threshold for missing data
thr <- 0.5
#names of the genes in subsets
APG <- c("PSMB8", "PSMB9", "TAP1", "TAP2", "TAPBP")
TAP <- c("TAP1", "TAP2")
PSMB <- c("PSMB8", "PSMB9")
APG_All <- c(APG, "APG", "PSMBs", "TAPs")
nonAPG <- c("BRD2", "DAXX", "KIFC1", "RGL2", "RXRBA")
APGmean <- c("APG")
nonAPGmean <- c("nonAPG")
MHC <- c("MHCI")
ex3 <- c("ex3")
g_colors <- c("PSMB8" = "darkorange", "PSMB9" = "darkorange4", "TAP1" = "green", "TAP2" = "green4", "TAPBP" = "mediumorchid1",
              "BRD2" = "gray60", "DAXX" = "turquoise1", "KIFC1" = "blue", "RGL2" = "gray30", "RXRBA" = "coral",
              "APG" = "blue", "nonAPG" = "red", "MHCI" = "magenta")

# FUNCTIONS #

#subsets sum_div to get genes of interest for a model
#prepares df which has separate columns for PDs used as 
#response, explanatory and "covariate" variables
#this is ready for PGLS modeling
div_to_xy <- function(arg_list){
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
  d <- sum_div %>% ungroup() %>% filter(gene_sym %in% xyz)
  #get the response as separate df
  dy <- d %>% filter(gene_sym %in% y)
  #seed the result with response
  dxyz <- dy
  for (i in xz) {
    #add gene name for PDs variables
    alf <- sym(paste0("PD_alpha_", i))
    gam <- sym(paste0("PD_gamma_", i))
    di <- d %>% filter(gene_sym %in% i) %>% select(sp_abr1, seq_type, q, PD_alpha, PD_gamma) %>% 
      rename(!!alf := PD_alpha, !!gam := PD_gamma)
    #put everything together ready for PGLSmodeling
    dxyz <- dxyz %>% left_join(di, by = c("sp_abr1", "seq_type", "q"))
  }
  return(dxyz)
}

#takes a list in the form (x = string, y = string, z(optional) = string)
#fits PGLS, returns df containing various outputs including the model results in themselves
#l_bounds is a numeric vecor (upper, lower), needed because sometimes the algorithm fails to converge
#then moving limits away from the boundaries sometimes helps
fit_gls_caper <- function (xyz, l_bounds, phylo = tree, transf = "") {
  d <- div_to_xy(xyz)
  x <- get(xyz[["x"]])
  y <- get(xyz[["y"]])
  z <- if (!is.null(xyz[["z"]])) get(xyz[["z"]]) else NULL
  #print(z)
  #not optimal but the list is not very long
  res <- NULL
  #iterates through dependent variables (commonly only one) 
  for (yi in y) {
    #iterates through sequence types, 
    #in this version DNA is left out
    for (st in c("codon", "AA", "AAGhm")) {
      #iterates through q values
      for (qval in c(0:2)){
        #get the right subset of observations
        di <- d %>% filter(gene_sym == yi, seq_type == st, q == qval) %>% as.data.frame()
        if (transf != "") {
          di <- di %>% mutate(PD_alpha = ifelse(PD_alpha >= 0.001, PD_alpha, 0.001),
                        PD_gamma = ifelse(PD_gamma >= 0.001, PD_gamma, 0.001))
        }
        #create comparative data object
        #uses the tree from the global environment
        di_comp <- comparative.data(phylo, di, taxon, vcv = TRUE)
        #iterates through alpha and gemma
        for (div_i in c("alpha", "gamma")) {
          #two explanatory variables
          if (!is.null(z)) {
            dep <- paste0(transf, "(PD_", div_i, ")")
            e <- paste0(transf, "(PD_", div_i, "_", x, ")")
            e_cov <- paste0(transf, "(PD_", div_i, "_", z, ")")
            f_int <- paste0(dep, " ~ ", e, " * ", e_cov)
            f_main <- paste0(dep, " ~ ", e, " + ", e_cov)
            print(c(yi, st, qval, div_i))
            tryCatch(m_int <- pgls(as.formula(f_int), data = di_comp, lambda = 'ML',  bounds=list(lambda=l_bounds)),
                      error = function(e){
                        message("An_error occurred:\n", e)
                        m_int <- pgls(as.formula(f_int), data = di_comp, lambda = 1)
                      })
            tryCatch(m_main <- pgls(as.formula(f_main), data = di_comp, lambda = 'ML', bounds=list(lambda=l_bounds)),
                      error = function(e){
                        message("An_error occurred:\n", e)
                        m_main <- pgls(as.formula(f_main), data = di_comp, lambda = 1)
                      })
            #m_main <- pgls(as.formula(f_main), data = di_comp, lambda = 'ML', bounds=list(lambda=l_bounds))
            dAICc <- m_int$aicc - m_main$aicc
            t_table <- as.data.frame(summary(m_main)$coefficients)
            p <- t_table$`Pr(>|t|)`
            est <- t_table$`Estimate`
            lam <- summary(m_main)$param[["lambda"]]
            mi_res <- data.frame(form = f_main, x = x, y = yi, cov = z, 
                                 seq_type = st, div_type = div_i, q = qval,
                                 Py = p[[2]], Pcov = p[[3]], dAICc_main_best = dAICc,
                                 y_est = est[[2]], cov_est = est[[3]], lambda = lam,
                                 mod = I(list(list(main = m_main, interaction = m_int))))
            
          #one explanatory variable
          } else {
            dep <- paste0(transf, "(PD_", div_i, ")")
            e <- paste0(transf, "(PD_", div_i, "_", x, ")")
            f <- paste0(dep, " ~ ", e)
            print(c(yi, st, qval, div_i))
            tryCatch(m_main <- pgls(as.formula(f), data = di_comp, lambda = 'ML', bounds=list(lambda=l_bounds)),
                      error = function(e){
                        message("An_error occurred:\n", e)
                        m_main <- pgls(as.formula(f), data = di_comp, lambda = 1)
                      })
            #m_main <- pgls(as.formula(f), data = di_comp, lambda = 'ML', bounds=list(lambda=l_bounds))
            t_table <- as.data.frame(summary(m_main)$coefficients)
            p <- t_table$`Pr(>|t|)`
            est <- t_table$`Estimate`
            lam <- summary(m_main)$param[["lambda"]]
            #currently only 
            mi_res <- data.frame(form = f, x = x, y = yi, cov = NA, seq_type = st, div_type = div_i, q = qval,
                                 Py = p[[2]], Pcov = NA, dAICc_main_best = NA,
                                 y_est = est[[2]], cov_est = NA, lambda = lam,
                                 mod = I(list(main = m_main)))
            
          }
          res <- rbind(res, mi_res)
        }
      }
    }
  }
  return(res)
}

#BODY#

#Read data and metadata ####

taxonomy <- read_tsv("taxonomy.txt")
genes <- read.table("gene_ids.txt", sep = "\t", header = T, stringsAsFactors = F)
cds <- readRDS("CDS_lengths.rds")
brd_mhc_al <- readRDS("Funct_alleles_df_BRD2_MHC_ex2_ex3_all_genera.rds")
tree <- read.nexus("30_salamanders_timetree.tre")

for (suff in c("", "_15_best_covered", "classical", "classical_and_intermediate")) {
  #"classical_and_intermediate" and "classical" use the full MIP dataset, 
  #they diverge only in MHC datasets
  if (suff == "classical_and_intermediate" || suff == "classical") {
    d <- readRDS(paste0("All_sp_div_excluding_stop_min20cov.rds"))
  } else {
    d <- readRDS(paste0("All_sp_div_excluding_stop_min20cov", suff, ".rds"))
  }
  m <- readRDS(paste0("All_sp_div_MHC_BRD", suff, ".rds"))
  div_sum_txt <- paste0("Diversities_summary_long_format", suff, ".txt")
  div_sum_rds <- paste0("Diversities_summary_long_format", suff, ".rds")
  pgls_file_name <- paste0("PGLS_modeling_results", suff)
  m <- m %>% rename(sp_abr1 = species) %>% left_join(taxonomy, by = "sp_abr1")

  #Per gene/group diversities ####

  #add gene symbols and cds lengths
  #consider remove unneccesary columns
  divAPGnAPG <- d %>% left_join(genes, by = c("gene")) %>% left_join(cds, by = c("species", "gene_sym"))

  #calculate diversities for each gene
  #these are sums over segments weighted by segment length

  sum_divAPGnAPG <- divAPGnAPG %>% filter(fr_typed >= thr) %>%
    group_by(species, class, gene_sym, seq_type, q) %>%
    mutate(seg_weight = len/sum(len),
           PD_alpha_weight = PD_alpha*seg_weight,
           PD_gamma_weight = PD_gamma*seg_weight) %>%
    summarise(reseq_len = sum(len),
              PD_alpha = sum(PD_alpha_weight),
              PD_gamma = sum(PD_gamma_weight)) %>%
    rename(sp_abr1 = species) %>% left_join(taxonomy, by = "sp_abr1")

  #arithmetic averages for PSMBs and TAPs
  psmb89 <- sum_divAPGnAPG %>% ungroup() %>% filter(gene_sym %in% PSMB) %>% select(-gene_sym) %>%
    group_by(sp_abr1, class, seq_type, q, family, genus, species, gen_abr1, sp_abr2) %>%
    summarise(reseq_len = sum(reseq_len),
              PD_alpha = mean(PD_alpha, na.rm = TRUE),
              PD_gamma = mean(PD_gamma, na.rm = TRUE)) %>%
    mutate(gene_sym = "PSMBs")

  tap12 <- sum_divAPGnAPG %>% ungroup() %>% filter(gene_sym %in% TAP) %>% select(-gene_sym) %>%
    group_by(sp_abr1, class, seq_type, q, family, genus, species, gen_abr1, sp_abr2) %>%
    summarise(reseq_len = sum(reseq_len),
              PD_alpha = mean(PD_alpha, na.rm = TRUE),
              PD_gamma = mean(PD_gamma, na.rm = TRUE)) %>%
    mutate(gene_sym = "TAPs")

  #combine the three above
  sum_divAPGnAPG <- bind_rows(sum_divAPGnAPG, psmb89, tap12)
  #unique(sum_divAPGnAPG[["gene_sym"]])

  #process BRD2 and MHC
  #for MHC involves weighted average across both exons

  #get mean allele length for BRD2 and MHC exons
  all_len <- brd_mhc_al %>% group_by(genus, species, exon) %>%
    summarise(mean_len_bp = mean(len, na.rm = TRUE))

  divBRD2MHC <- m %>% left_join(all_len, by = c("genus", "species","exon"))
  to_del <- c("exon", "n_hap", "N_typed", "S", "S_unamb", "dmax", "PD_beta", "local_similarity", "region_similarity")

  divBRD2 <- divBRD2MHC %>% filter(exon == "brd2") %>% select(-all_of(to_del)) %>%
    add_column(class = "nonAPG", gene_sym = "BRD2", .after = "sp_abr1") %>%
    relocate(mean_len_bp, PD_alpha, .after = "q") %>%
    rename(reseq_len = mean_len_bp)

  ##check whether names match
  #identical(sort(names(sum_divAPGnAPG)), sort(names(divBRD2)))

  divMHC <- divBRD2MHC %>% filter(exon != "brd2") %>%
    select(-all_of(to_del))  %>%
    group_by(sp_abr1, family, genus, species, gen_abr1, sp_abr2, seq_type, q) %>%
    mutate(seg_weight = mean_len_bp/sum(mean_len_bp),
           PD_alpha_weight = PD_alpha*seg_weight,
           PD_gamma_weight = PD_gamma*seg_weight) %>%
    summarise(reseq_len = sum(mean_len_bp),
              PD_alpha = sum(PD_alpha_weight),
              PD_gamma = sum(PD_gamma_weight)) %>%
    mutate(gene_sym = "MHCI",
           class = "MHC")
  
  divex3 <- divBRD2MHC %>% filter(exon == "3") %>%
    select(-all_of(to_del))  %>%
    group_by(sp_abr1, family, genus, species, gen_abr1, sp_abr2, seq_type, q) %>%
    mutate(seg_weight = 1,
           PD_alpha_weight = PD_alpha*seg_weight,
           PD_gamma_weight = PD_gamma*seg_weight) %>%
    summarise(reseq_len = mean_len_bp,
              PD_alpha = sum(PD_alpha_weight),
              PD_gamma = sum(PD_gamma_weight)) %>%
    mutate(gene_sym = "ex3",
           class = "MHC")

  sum_div <- bind_rows(sum_divAPGnAPG, divMHC, divBRD2, divex3) %>% select(-c(gen_abr1, sp_abr2)) %>% ungroup()

  #calculates average diversities for APGs and nonAPGs
  temp <- sum_div %>% group_by(sp_abr1, class, seq_type, q, family, genus, species) %>%
    summarise(reseq_len = sum(reseq_len),
              PD_alpha = mean(PD_alpha, na.rm = TRUE),
              PD_gamma = mean(PD_gamma, na.rm = TRUE)) %>%
    filter(class != "MHC") %>%
    mutate(gene_sym = class)

  #adds average diversities calculated above to the diversity table
  #this is the complete diversity table
  #can be extended if needed, e.g. with weighted averages
  sum_div <- bind_rows(sum_div, temp) %>% mutate(taxon = paste0(genus, "_", species))

  write_lnx_head(sum_div, div_sum_txt)
  saveRDS(sum_div, div_sum_rds)

  #PGLS modeling ####
   
  #list of variables for modeling with fit_gls_caper
  al <- list(list(x = "MHC", y = "APG_All", z = "nonAPGmean"),
               list(x = "ex3", y = "APG_All", z = "nonAPGmean"),
               list(x = "nonAPGmean", y = "MHC"),
               list(x = "nonAPGmean", y = "ex3"))
  #l_bound may need to be adjusted if the algorithm failed to converge 
  l <- lapply(al, fit_gls_caper, l_bounds = c(0.001, 0.999), phylo = tree)
  l_all <- bind_rows(l)
  saveRDS(l_all, paste0(pgls_file_name, ".rds"))
  write_lnx_head(select(l_all, -mod), paste0(pgls_file_name, ".txt"))
}




# #Various tests and checks ####
# 
# sum_div_flt <- sum_div %>% filter(gene_sym %in% c("PSMBs", "TAPs", "MHCI", "APG", "nonAPG"))
# sum_div_flt <- sum_div %>% filter(gene_sym %in% c("APG", "nonAPG"))
# #plot distributions of variables
# ggplot(sum_div_flt, aes(x = PD_alpha, fill = as.factor(q))) +
#    geom_histogram(bins = 7, position = "dodge") +
#    facet_grid(rows = vars(seq_type), cols = vars(gene_sym), scales = "free")
# 
# 
# # #plot distributions of variables
# ggplot(sum_div_flt, aes(sample = PD_gamma, colour = as.factor(q))) + 
#   stat_qq() + stat_qq_line() +
#   facet_grid(rows = vars(seq_type), cols = vars(gene_sym), scales = "free")
# #ggsave("alpha.pdf")
# #   geom_histogram(bins = 7, position = "dodge") + 
# 
# # #plot distributions of variables
# ggplot(sum_div_flt, aes(sample = log(PD_gamma), colour = as.factor(q))) + 
#   stat_qq() + stat_qq_line() +
#   facet_grid(rows = vars(seq_type), cols = vars(gene_sym), scales = "free")
# #ggsave("alpha.pdf")
# 



#Make plots ####

sum_div <- readRDS("Diversities_summary_long_format.rds") %>% filter(seq_type != "DNA") %>% 
  mutate(seq_type = ifelse(seq_type == "codon", "dS", seq_type))

arg_to_plot <- list(list(x = "MHC", y = "APG"),
                    list(x = "MHC", y = "nonAPG"),
                    list(x = "nonAPGmean", y = "APG"),
                    list(x = "nonAPGmean", y = "MHC"))

plot_div <- function(arg_list){
  d <- div_to_xy(arg_list)
  x <- get(arg_list[["x"]])
  adj_plot <- function(p){
    adj_p <- p + geom_point() + geom_smooth(method='lm', se = FALSE) +
      scale_colour_manual(values = g_colors) +
      facet_grid(rows = vars(seq_type), cols = vars(q), scales = "free") +
      theme_bw()
    print(adj_p)

  }
  a_plot <- ggplot(d, aes(x = !!sym(paste0("PD_alpha_", x)), y = PD_alpha, color = gene_sym)) %>% adj_plot()
  g_plot <- ggplot(d, aes(x = !!sym(paste0("PD_gamma_", x)), y = PD_gamma, color = gene_sym)) %>% adj_plot()
  a_g <- ggarrange(a_plot, g_plot,
                         labels = c("alpha", "gamma"),
                         ncol = 2, nrow = 1,
                         font.label = list(size = 11, color = "black", face = "bold", family = NULL))
  ggsave(paste0(arg_list[["y"]], "_vs_", arg_list[["x"]], ".pdf"), width = 20, height = 10, units = "in")
}
lapply(arg_to_plot, plot_div)

