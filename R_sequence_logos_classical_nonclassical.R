library(ggseqlogo)
library(tidyverse)
library(ggpubr)
library(ape)
library(phangorn)
library(treeio)
library(ggtree)
library(seqRFLP)

source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

# FUNCTIONS ####

#takes list of positions in aa alignment
#with expected aa (accepts regex)
#counts the number of matching positions
#the idea was to count the number of expected "anchor residues"
nclass <- function(s, pos, aa) {
  ss <- str_split(s, "") %>% unlist() %>% .[pos]
  score <- 0
  for(e in seq_along(ss)) {
    m <- grepl(aa[e], ss[e])
    score <- score + m
  }
  return(score)
}
#vectorize but only the string with sequence, two other args are nit vectorizes
#USE.NAMES = FALSE, because otherwise there's a problem with returning named vector, not sure why
vect_nclass <- Vectorize(nclass, vectorize.args = c("s"), USE.NAMES = FALSE)

#adaptation of genus_trees() that doesn't mark species but marks 
#potential classical/nonclassical provided as a column of gen
genus_trees_class_nonclass <- function(t, genotypes = gen){
  d <- genotypes %>% filter(genus == t)
  #c_n <- d %>% pull(non_classical) %>% unique() %>% sort()
  #sp_list <- d %>% pull(species) %>% unique() %>% sort()
  ts <- vector("list", 2)
  for(ex in c("2", "3")){
    seq <- d %>% filter(exon == ex) %>% select(label, seq) %>% distinct(label, seq) %>% as.data.frame() %>% df2DNA()
    ann <- d %>% filter(exon == ex) %>% select(label, funct, genus, species, class_nonclass) %>% 
      distinct(label, species, .keep_all = TRUE) %>% mutate(pres = 1) %>%
      pivot_wider(names_from = "species", values_from = "pres", values_fill = list(pres = 0))
    alg_file <- paste0(t, "_ex_", ex, "_alignment.rds")
    if(file.exists(alg_file)){
      alg <- readRDS(alg_file)
    } else {
      alg <- muscle(seq, exec = "muscle3.8.31_i86win32.exe")
      saveRDS(alg, alg_file)
    }
    JC <- dist.dna(alg, model = "JC", pairwise.deletion = TRUE)
    tree <- bionj(JC)
    tree <- midpoint(tree)
    t_ann <- full_join(tree, ann, by = "label")
    t_plot <- ggtree(t_ann, layout = "rectangular", ladderize = TRUE, aes(color = funct), size = 0.5) + 
      scale_colour_manual(name = NULL, values = c("y" = "green", "n" = "red"), 
                          labels = c("y" = "funct", "n" = "non-funct"),
                          na.value = "darkgrey") +
      geom_tippoint(aes(fill = class_nonclass), size = 1.5, shape = 21, colour = "white") +
      scale_fill_manual(name = element_blank(),
                        values = c("class" = "green", "intermed" = "yellow", "nonclass" = "red"),
                        labels = c("class" = "classical", "intermed" = "intermediate", "nonclass" = "non-classical"),
                        na.value = "lightgrey")
    theme_tree2(legend.title = element_blank()) 
    ts[[as.character(ex)]] <- t_plot
  }
  
  tree_ex23 <- ggarrange(ts[["2"]], ts[["3"]],
                         labels = c("ex2", "ex3"),
                         ncol = 2, nrow = 1,
                         font.label = list(size = 11, color = "black", face = "bold", family = NULL),
                         common.legend = TRUE)
  
  Fig <- annotate_figure(tree_ex23,
                         top = text_grob(paste0(t), face = "italic", size = 14))
  return(Fig)
}

# BODY

#Generate sequence logo ####

taxonomy <- read_tsv("taxonomy.txt") %>% select(family, genus, gen_abr1)

ex2aa <- read.FASTA("sequences/All_genera_ex2_funct_full_codons.fas") %>% trans() %>% bin2df()
ex3aa <- read.FASTA("sequences/All_genera_ex3_funct_full_codons.fas") %>% trans() %>% bin2df()
ex2_ex3_aa_df <-bind_rows(ex2aa, ex3aa)  %>% rename(seq_aa = seq)
saveRDS(ex2_ex3_aa_df, "ex2_ex3_aa_alignment_df.rds")

s <- ex2_ex3_aa_df %>% 
  separate(label, into = c("gen_abr1", "exon", "allele"), sep = "_") %>% 
  left_join(taxonomy, by = "gen_abr1") %>% 
  mutate(fam_gen = paste0(family, ": ", genus))

genera <- s %>% pull(fam_gen) %>% unique() %>% sort()

lex2 <- list()
lex3 <- list()
for (g in genera) {
  l <- s %>% filter(fam_gen == g, exon == "ex2") %>% pull(seq_aa)
  lex2[[g]] <- l
  l <- s %>% filter(fam_gen == g, exon == "ex3") %>% pull(seq_aa)
  lex3[[g]] <- l
}

#anchor positions in ALIGNMENT of each exon
ex2_anchor <- c(61)
ex3_anchor <- c(47, 50, 51, 66)

#PBSes common to all HLA-A-C
ex2_pbs_all_pol <- c(2, 27, 65, 68, 69, 72, 75, 76, 79, 82, 83)
ex3_pbs_all_pol <- c(14, 23, 24, 58, 62, 70)
#PBS in some HLA
ex2_pbs_some_pol <- c(4, 16, 18, 40)
ex3_pbs_some_pol <- c(11, 46, 71)

mytheme <- theme(axis.text.x = element_text(angle = 90, size = 4.5),
      strip.text = element_text(size = 8, colour = "gray30"),
      plot.title = element_text(hjust = 0.5, size = 14),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "none",
      panel.spacing = unit(0, "lines"),
      strip.text.x = element_text(margin = margin(0,0,0,0, "cm")))

ex2_logo <- ggseqlogo(lex2, ncol = 1, method = 'prob', facet = "wrap", scales = "fixed") + 
  ggtitle("exon 2") +
  annotate('rect', xmin = ex2_anchor - 0.5, xmax = ex2_anchor + 0.5, ymin = -0.2, ymax = 1.2, 
           alpha = .12, col=NA, fill='blue') +
  annotate('rect', xmin = ex2_pbs_all_pol - 0.5, xmax = ex2_pbs_all_pol + 0.5, ymin = -0.2, ymax = 1.2, 
           alpha = .4, col=NA, fill='gold') +
  annotate('rect', xmin = ex2_pbs_some_pol - 0.5, xmax = ex2_pbs_some_pol + 0.5, ymin = -0.2, ymax = 1.2, 
           alpha = .2, col=NA, fill='yellow') +
  mytheme

ex3_logo <- ggseqlogo(lex3, ncol = 1, method = 'prob', facet = "wrap", scales = "fixed") + 
  ggtitle("exon 3") +
  annotate('rect', xmin = ex3_anchor - 0.5, xmax = ex3_anchor + 0.5, ymin = -0.2, ymax = 1.2, 
           alpha = .12, col=NA, fill='blue') +
  annotate('rect', xmin = ex3_pbs_all_pol - 0.5, xmax = ex3_pbs_all_pol + 0.5, ymin = -0.2, ymax = 1.2, 
           alpha = .4, col=NA, fill='gold') +
  annotate('rect', xmin = ex3_pbs_some_pol - 0.5, xmax = ex3_pbs_some_pol + 0.5, ymin = -0.2, ymax = 1.2, 
           alpha = .2, col=NA, fill='yellow') +
  mytheme

ex2_ex3_logo <- ggarrange(ex2_logo, ex3_logo, ncol = 2, widths = c(1.05, 1))
ggsave("Fig_S1_ex2_ex3_logo.pdf", width = 28, height = 15, units = "cm")

#Identification of putative classical, intermediate and nonclassical alleles ####
#Alleles with the canonical aa in all anchor residues formed "conserved_anchor dataset

#positions in alignment with anchor residues
pos_ex2 <- c(61)
pos_ex3 <- c(47, 50, 51, 66)
#residues in these positions
aa_ex2 <- c("Y")
aa_ex3 <- c("T", "[KR]", "W", "Y")

#reads data
d <- ex2_ex3_aa_df %>% 
  separate(label, into = c("gen_abr", NA, NA), sep ="_", remove = FALSE)
ex2 <- d %>% filter(grepl("_ex2_", label)) %>% 
  mutate(n_class_aa = vect_nclass(seq_aa, pos_ex2, aa_ex2),
         class_nonclass = ifelse(n_class_aa == 0, "nonclass", "class"))
ex3 <- d %>% filter(grepl("_ex3_", label)) %>% 
  mutate(n_class_aa = vect_nclass(seq_aa, pos_ex3, aa_ex3),
         class_nonclass = ifelse(n_class_aa == 4, "class", ifelse(n_class_aa ==3, "intermed", "nonclass")))
ex2_ex3 <- bind_rows(ex2, ex3)
#saveRDS(select(ex2_ex3, label, class_nonclass), "ex2_ex3_class_nonclass_status.rds")


# #summary of numbers and fractions of non-classical
# ex2_sum <- ex2 %>% group_by(gen_abr) %>%
#   summarise(n_class = sum(n_class_aa == 1),
#             fr_class = n_class/n())
# ex3_sum <- ex3 %>% group_by(gen_abr) %>%
#   summarise(n_class_min3aa = sum(n_class_aa >= 3),
#             n_class_4aa = sum(n_class_aa == 4),
#             fr_class_min3aa = n_class_min3aa/n(),
#             fr_class_4aa = n_class_4aa/n())
# ex2_sum_sum <- ex2 %>%
#   summarise(n_class = sum(n_class_aa == 1),
#             fr_class = n_class/n())
# ex3_sum_sum <- ex3 %>%
#   summarise(n_class_min3aa = sum(n_class_aa >= 3),
#             n_class_4aa = sum(n_class_aa == 4),
#             fr_class_min3aa = n_class_min3aa/n(),
#             fr_class_4aa = n_class_4aa/n())

# #reads genotypes and adds info about classicity
# g <- readRDS("genotypes_ex2_ex3_brd2_all_genera.rds")
# ex2_ex3_annot <- ex2_ex3 %>% select(label, class_nonclass)
# gen <- g %>% left_join(ex2_ex3_annot, by = "label")
# 
# #generate trees for genera
# genera <- gen %>% pull(genus) %>% unique()
# genera_trees <- lapply(genera, genus_trees_class_nonclass)
# 
# #plot trees for genera to a single pdf
# pdf(file = "Trees_ex2_ex3_genera_classical_nonclassical.pdf", width = 10, height = 8, onefile = TRUE)
# for(i in seq_along(genera_trees)){
#   plot(genera_trees[[i]])
# }
# dev.off()

# #generate alignments of only classical and classical plus intermediate alleles
# #for exon 3, for all species and each genus separately
# alg <- readRDS("Funct_alleles_df_BRD2_MHC_ex2_ex3_all_genera.rds")
# status <- readRDS("ex2_ex3_class_nonclass_status.rds")
# alg_ex3 <- alg %>% filter(exon == "3") %>% left_join(status, by = "label")
# ex3_class <- alg_ex3 %>% filter(class_nonclass == "class") %>% select(genus, label, seq_codon)
# ex3_class_intermed <- alg_ex3 %>% filter(class_nonclass == "class" | class_nonclass == "intermed") %>% 
#   select(genus, label, seq_codon)
# cl <- ex3_class %>% select(label, seq_codon) %>% as.data.frame()
# dataframe2fas(cl, "Ex_3_classical_all_genera.fas")
# cl_int <- ex3_class_intermed %>% select(label, seq_codon) %>% as.data.frame()
# dataframe2fas(cl_int, "All_genera_ex_3_classical_and_intermediate.fas")
# genera <- ex3_class %>% pull(genus) %>% unique()
# 
# for (g in genera){
#   cl <- ex3_class %>% filter(genus == g) %>% select(-genus) %>% as.data.frame()
#   cl_int <- ex3_class_intermed %>% filter(genus == g) %>% select(-genus) %>% as.data.frame()
#   cl_name <- paste0(g, "_ex3_classical.fas")
#   cl_int_name <- paste0(g, "_ex3_classical_and_intermediate.fas")
#   dataframe2fas(cl, cl_name)
#   dataframe2fas(cl_int, cl_int_name)
# }
# 