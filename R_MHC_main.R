#load packages and functions ####
library(ape)
library(treeio)
library(ggtree)
library(phangorn)
library(ggpubr)
library(hillR)
library(seqRFLP)
library(phytools)
library(tidyverse)


source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

# F U N C T I O N S #

#reads genotypes for a genus from text files
#MHC exons and BRD
read_gen <- function(tax){
  print(tax)
  gen_tax <- paste0("genotypes/", tax)
  ex2 <- read_tsv(paste0(gen_tax, "_ex2.txt")) %>% 
    pivot_longer(-seq, names_to = "id", values_to = "n_reads") %>%
    mutate(exon = "2")
  ex3 <- read_tsv(paste0(gen_tax, "_ex3.txt")) %>% 
    pivot_longer(-seq, names_to = "id", values_to = "n_reads") %>%
    mutate(exon = "3")
  brd <- read_tsv(paste0(gen_tax, "_BRD2_genotypes.txt")) %>%
    pivot_longer(-seq, names_to = "id", values_to = "n_reads") %>%
    mutate(exon = "brd2")
  ex2_3_brd2 <- bind_rows(ex2, ex3, brd) %>% 
    filter(n_reads > 0) %>% mutate(id = as.numeric(id)) %>% select(id, exon, seq, n_reads)
  return(ex2_3_brd2)
}

#takes dataframe (df) with genotypes in long format
#and exon id (e)
#returns list with named elements:
#"DNA", "codon", "AA" - each of these is a list of three elements:
#1) df with binary encoded genotypes (inds in rows, alleles in cols), col1: individual id
#2) a single row df with basic statistics of the segment
#3) dataframe with allele sequences
#4) dataframe of allele incidences - for each allele number of individuals in which it was present
#IMPORTANT! renames alleles according to abundance
bin_gen_MHCBRD <- function(e, df){
  sp <- df %>% pull(sp_abr1) %>% unique()
  dff <- df %>% filter(exon == e) %>% select(id, label, seq, seq_codon, seq_aa, presence)
  res <- NULL
  for(s in c("codon", "AA")){
    if(s == "codon") h <- "seq_codon"
    if(s == "AA") h <- "seq_aa"
    a <- dff %>% group_by(across(h)) %>% summarize(n_ind = n()) %>% 
      arrange(desc(n_ind), desc(h)) %>% mutate(allele = paste0("all_", sprintf("%02d",row_number())),
                                               presence = 1)
    seq <- a %>% select(allele, h) %>% arrange(allele) %>% as.data.frame()
    if(s == "AA"){
      seq <- df2prot(seq)
    } else {
      seq <- df2DNA(seq)
    }
    #allele incdence - number of individuals having the allele
    #useful for drawing "allele frequency spectra"
    inc <- a %>% select(allele, n_ind) %>% 
      mutate(species = sp,
             exon = e) %>% select(species, exon, allele, n_ind) %>% 
      as.data.frame()
    d <- dff %>% select(id, h)
    print(d)
    d_a <- d %>% left_join(a, by = h) %>% select(-c(n_ind, h)) %>% distinct() %>% 
      arrange(allele) %>% 
      pivot_wider(names_from = allele, values_from = presence, values_fill = list(presence = 0)) %>%
      arrange(id) %>% as.data.frame()
    print(d_a)
    n_hap <- ncol(d_a) - 1
    n_ind <- nrow(d_a)
    summary <- data.frame("species" = sp, "exon" = e, "n_hap" = n_hap, 
                          "N_typed" = n_ind)
    rownames(summary) <- c()
    r <- list("genotypes" = d_a, "summary" = summary, "seq" = seq, "all_incidence" = inc)
    print(d_a$id)
    res[[s]] <- r
  }
  res[["AAGhm"]] <- res[["AA"]]
  res[["DNA"]] <- res[["codon"]]
  return(res)
}


#calculates alpha and gamma diversities for a species
#and alpha diversities for all individuals
#uses calc_div defined in R_functions_WB_MHC.R
get_div_MHCBRD <- function(t, fr = gen, alleles = fun_aln){
  #takes functional alleles from the focal species  
  t_df <- fr %>% select(sp_abr1, id, exon, label, funct, seq, seq_codon, seq_aa) %>% 
    filter(sp_abr1 == t & funct == "y") %>% 
    mutate(presence = 1)
  #get names of "exons"
  exons <- t_df %>% pull(exon) %>% unique()
  #to keep results of calc_div_MHCBRD
  dd <- vector("list", 3)
  #to keep results of bin_gen_MHCBRD
  bgs <- NULL
  #loop through exons
  #for each generate binary genotypes and calculate diversities
  #results in two lists initialised above
  Ghm <- read.csv("Grantham_Distance_table.csv", header = T, row.names = 1, encoding = "UTF-8")
  for(i in seq_along(exons)){
    exon <- exons[[i]]
    bg <- bin_gen_MHCBRD(exon, t_df)
    bgs[[exon]] <- bg
    #"DNA" is not included here in calculations, although it could be
    for(z in c("codon", "AA", "AAGhm")){
      d <- calc_div(bg[[z]], t, typ = z, Ghm)
      dd[[i]][[z]] <- d
    }
  }
  alpha_gamma <- NULL
  individual_alpha <- NULL
  #collect alpha_gamma and individual_alpha for all exons into lists
  for(z in c("DNA", "codon", "AA", "AAGhm")){
    seq_type <- data.frame(seq_type = z)
    #list, for each exon results for the type
    b <- lapply(dd, "[[", z)
    #list, for each exon species level alpha & gamma
    a_g <- lapply(b, "[[", "alpha_gamma")
    a_g <- bind_rows(a_g)
    #list, for each exon individual alphas
    i_a <- lapply(b, "[[", "individual_alpha")
    i_a <- bind_rows(i_a)
    alpha_gamma[[z]] <- bind_cols(seq_type, a_g)
    individual_alpha[[z]] <- bind_cols(seq_type, i_a)
  }
  #bind lists into df
  sp_div <- bind_rows(alpha_gamma)
  ind_div <- bind_rows(individual_alpha)
  saveRDS(sp_div, paste0("rds/", t, "_sp_div_MHC_BRD.rds"))
  saveRDS(ind_div, paste0("rds/", t, "_indivdual_alphadiv_MHC_BRD.rds"))
  saveRDS(bgs, paste0("rds/", t, "_binary_genotypes_MHC_BRD.rds"))
  res = list("sp_div" = sp_div, "ind_div" = ind_div, "bgs" = bgs)
  return(res)
} 

#constructs ex2 & 3 trees for a taxon, color codes functional/nonfunctional
#and species
#returns ready to plot object
genus_trees <- function(t, genotypes = gen){
  d <- genotypes %>% filter(genus == t)
  sp_list <- d %>% pull(species) %>% unique() %>% sort()
  ts <- vector("list", 2)
  for(ex in c("2", "3")){
    seq <- d %>% filter(exon == ex) %>% select(label, seq) %>% distinct(label, seq) %>% as.data.frame() %>% df2DNA()
    ann <- d %>% filter(exon == ex) %>% select(label, funct, genus, species) %>% 
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
      theme_tree2(legend.title = element_blank()) 
    col <- c("red", "blue", "chartreuse3", "cyan", "darkgoldenrod1", "darkorchid", "darkorange", "grey50", "pink")
    for(s in seq_along(sp_list)){
      print(s)
      to_add <- paste0("geom_tippoint(aes(subset = ", sp_list[[s]], " == 1, x = x + (0.02*", s, ")*max(x), fill = \"", 
                       sp_list[[s]], "\"), size = 1.5, shape = 21, color = \"white\")")
      to_add_fill <- paste0("scale_fill_manual(name = \"species\",
                          values = col,
                          guide = guide_legend(label.theme = element_text(angle = 0, face = \"italic\")))")
      print(to_add)
      t_plot <-  t_plot + eval(parse(text = to_add)) + eval(parse(text = to_add_fill))
    }
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

#identifies putatively nonfunctional MHC alleles within species
#these are all nonfunctional alleles designated on the basis of frameshifts/stop codons
#plus all alleles falling into the most inclusive clade (threshold given by thr <50; 100>)
#with at least one such allele
#returns vector with labels of all potentially nonfunctional alleles
#tax: genus
#nfun: labels and nonfunctional alleles (frameshift and stops)
#gen: genotypes used by most functions in this script
#CAUTION - ALTHOUGH THIS FUNCTION WORKS CORRECTLY, IN SOME CASES IT RETURNS ALMOST EVEYTHING
#THIS IS RELATED TO PECULIAR TOPOLOGIES IN SOME GENERA, SUCH AS ICHTHYOSAURA
#Visual inspection showed as unreliable at threshold 90:
#Bat_ex2 - no alleles should be added to nonfunctional
#Eur_ex2 - only one allele should be added
#Ich_ex2 - has to be done manually
#Ple_ex3 - no alleles should be added to nonfunctional
#Pro_ex2 - has to be done manually
#to facilitate this the function produces also trees with nonfunctional (stop/frameshift highlighted)
potentially_nf <- function(t, thr) {
  nfun_cl <- function(nf, tbp, t){
    #iterate through nonfunctional alleles
    if (length(nf) == 0) tips <- NA else {
      nodes <- integer(length(nf))
      for (i in seq_along(nf)) {
        n <- nf[[i]]
        p <- tbp %>% filter(label == n) %>% pull(parent)
        bp <- 0
        #go up the tree until you find clade with support >= threshold or you reach "root"
        while (p != n && bp < thr) {
          r <- tbp %>% filter(node == p)
          n <- r$node
          p <- r$parent
          if (p != n) bp <- ifelse(!is.na(r$clade_bp), r$clade_bp, 0)
        }
        #if bp < threshold there's no node to include
        if (bp >= thr) nodes[[i]] <- n else nodes[[i]] <- NA
      }
      nodes <- unique(nodes[!is.na(nodes)])
      tips <- NULL
      print(nodes)
      for (node in nodes) {
        l <- extract.clade(t, node)$tip.label
        tips <- union(tips, l)
      }
      tips <- union(tips, nf)
      print(tips)
    }
    return(tips[!is.na(tips)])
  }
  d <- gen %>% filter(genus == t)
  pot_nf <- NULL
  N_nf <- list("2" = 0, "3" = 0)
  for (ex in c("2", "3")) {
    seq <- d %>% filter(exon == ex) %>% select(label, seq) %>% distinct(label, seq) %>% as.data.frame() %>% df2DNA()
    all_tax <- names(seq)
    nfun_tax <- intersect(nfun, all_tax)
    alg_file <- paste0(t, "_ex_", ex, "_alignment.rds")
    if(file.exists(alg_file)){
      alg <- readRDS(alg_file)
    } else {
      alg <- muscle(seq, exec = "muscle3.8.31_i86win32.exe")
      saveRDS(alg, alg_file)
    }
    #function to construct bionj
    JCbionj <- function(x) bionj(dist.dna(x, model = "JC", pairwise.deletion = TRUE))
    #midpoint root - important
    tree <- midpoint.root(JCbionj(alg))
    #get bootstrapped trees
    bstrees <- boot.phylo(tree, alg, JCbionj, trees = TRUE)$trees
    ## get support for each clade:
    #get tibble with correct bootstrap supports for clades
    #taking simple bootstrap values doesn't work because of the issues at the root
    clad <- prop.clades(tree, bstrees, rooted = FALSE)
    bp2 <- tibble(node=1:Nnode(tree) + Ntip(tree), clade_bp = clad)
    #this is input to function that picks all labels
    #that are in inclusive clades with nonfunctional alleles
    #specifically, this tibble allows to identify nodes with descendants to extract
    tree_bp <- full_join(as_tibble(tree), bp2, by="node")
    pdf(file = paste0(t, "_ex_", ex, "_bp_labels.pdf"), width = 10, height = 8, onefile = TRUE)
    layout(1)
    par(mar = rep(2, 4), cex = 0.6)
    plot(tree, font = 1)
    # drawSupportOnEdges(clad)
    nodelabels(clad)
    i <- which(tree$tip.label %in% nfun_tax)
    if (length(i) > 0) tiplabels(tree$tip.label[i], i, adj = 0, cex = 0.6)
    # edgelabels()
    dev.off()
    nf_ex <- nfun_cl(nfun_tax, tree_bp, tree)
    pot_nf <- union(nf_ex, pot_nf)
    N_nf[[ex]] <- length(nf_ex)
  }
  print(pot_nf)
  N <- data.frame(genus = t, N_nf_ex2 = N_nf[["2"]], N_nf_ex3 = N_nf[["3"]])
  return(list(pot_nf = pot_nf, N = N))
}


# B O D Y #

#Loads data and extracts various info ####

#get allele info
al <- read_tsv("allele_info.txt")

#get classical/nonclassical status
class_nonclass <- readRDS("ex2_ex3_class_nonclass_status.rds")

#get ids of potentially functional and nonfunctional alleles
nfun <- al %>% filter(funct == "n") %>% pull(label)
fun <- al %>% filter(funct == "y") %>% pull(label)

#get coverage
#this is available only for ex2 and ex3, for brd will be calculated
cov <- read_tsv("genotypes/coverage.txt", col_types = "ici")

#count na.id????
id <- read_tsv("individual_info.txt")

#get various taxonomic info
taxonomy <- read_tsv("taxonomy.txt")

#genera abbreviations
tax <- taxonomy %>% pull(gen_abr1) %>% unique()

#vector of families & genera
fam_gen <- taxonomy %>% select(family, genus) %>% distinct()

#species abbreviations
sp_abr1 <- taxonomy %>% pull(sp_abr1) %>% unique() 

#vector of genera
genera <- taxonomy %>% pull(genus) %>% unique()

#vector of families
families <- taxonomy %>% pull(family) %>% unique()

#read codon-based alignments of functional alleles
#these include everything except alleles with frameshifts/stop codons
ex2_fun_aln <- read.FASTA("sequences/All_genera_ex2_funct_full_codons.fas")
ex3_fun_aln <- read.FASTA("sequences/All_genera_ex3_funct_full_codons.fas")
brd2_fun_aln <- read.FASTA("sequences/All_genera_brd2_funct_full_codons.fas")
fun_aln <- list("2" = ex2_fun_aln,
                "3" = ex3_fun_aln,
                "brd2" = brd2_fun_aln)
fun_aa <- lapply(fun_aln, trans)
fun_codon_df <- lapply(fun_aln, bin2df) %>% bind_rows() %>% rename(seq_codon = seq)
fun_aa_df <- lapply(fun_aa, bin2df) %>% bind_rows() %>% rename(seq_aa = seq)
fun_codon_aa_df <- fun_codon_df %>% left_join(fun_aa_df, by = "label") 

#Reads or generates genotypes ####
gen_file <- "genotypes_ex2_ex3_brd2_all_genera.rds"
if(file.exists(gen_file)){
  print("I'm reading a pre-existing genotype .rds file!")
  gen <- readRDS(gen_file)
} else {
  #read genotypes for all genera into a list
  g <- lapply(tax, read_gen)
  taxonomy_sel <- taxonomy %>% select(family, genus, gen_abr1, species, sp_abr1)
  gen <- bind_rows(g) %>% 
    left_join(al, by = "seq") %>% 
    left_join(fun_codon_aa_df, by = "label") %>% 
    left_join(id, by = "id") %>% 
    left_join(taxonomy_sel, by = c("genus", "species"))
  brd_cov <- gen %>% filter(exon == "brd2") %>% group_by(id) %>% mutate(cov = sum(n_reads)) %>%
     select(id, exon, cov) %>% distinct()
  cov <- bind_rows(cov, brd_cov)
  gen <- gen %>% left_join(cov, by = c("id", "exon"))
  saveRDS(gen, gen_file)
}


#List of typed individuals ####
#ids of individuals genotyped in ex2, ex3 and BRD2
ids_gen <- gen %>% select(id, exon, sp_abr1) %>% distinct()
exons <- ids_gen %>% pull(exon) %>% unique()
typed <- NULL
for (s in sp_abr1) {
  temp1 <- ids_gen %>% filter(sp_abr1 == s)
  typed_ex <- NULL
  for (e in exons) {
    ids <- temp1 %>% filter(exon == e) %>% pull(id) %>% unique() %>% sort()
    typed_ex[[e]] <- ids
  }
  typed[[s]] <- typed_ex
}
saveRDS(typed, "Ids_typed_amplicons.rds")

# #Phylogeny-based potentially nonfunctional alleles ####
# #Getting list of all potentially nonfunctional alleles
# #including both frameshifts/stops and alleles present 
# #together with them in clades withmin 90% bootstrap support 
# #partly automated, partly manual
# #p <- lapply(genera, potentially_nf, thr = 90)
# #saveRDS(p, "potentially_nonfunctional_MHC.rds")
# ##get dataframe with counts of potentially nofunctional alleles
# ##useful for troubleshooting and identification of problematic genera
# #nf <- bind_rows(lapply(p, "[[", 2))
# 
# #get saved results of automated search
# temp <- readRDS("potentially_nonfunctional_MHC.rds") %>% 
#   lapply("[[", 1) %>% unlist()
# 
# #remove taxa where automatic assignment was incorrect
# to_rem <- c("Bat_ex2", "Eur_ex2", "Ich_ex2", "Kar_ex2", "Ple_ex3", "Pro_ex2")
# temp1 <- grep(paste(to_rem, collapse = "|"), temp, value = TRUE, invert = TRUE)
# 
# #for genera where automated protocol didn't work I extracted additional potentially nonfunctional alleles
# nf_Eur_Ich_Pro <- c("Eur_ex2_023", "Ich_ex2_021", "Pro_ex2_116", "Pro_ex2_010", "Pro_ex2_014",
#                     "Pro_ex2_020", "Pro_ex2_021", "Pro_ex2_052", "Pro_ex2_035", "Pro_ex2_097",
#                     "Pro_ex2_059", "Pro_ex2_030", "Pro_ex2_046", "Pro_ex2_058", "Pro_ex2_082", "Pro_ex2_018",
#                     "Pro_ex2_121")
# 
# #add_everything and make unique
# #these are nonfunctional and potentially nonfunctional alleles
# nf <- union(temp1, nf_Eur_Ich_Pro)
# 
# #and these are to be removed from alignments of potentially functional
# remove_from_fun <- intersect(nf, fun)
# saveRDS(remove_from_fun, "Alleles_to_remove_from_functional.rds")


#Get sample sizes of 15 ####

typed <- readRDS("Ids_typed_MIPs_amplicons.rds")
mhc15 <- lapply(typed, "[[", "mhc15") %>% unlist()
brd15 <- lapply(typed, "[[", "brd15") %>% unlist()
gen15 <- gen %>% filter((exon %in% c("2", "3") & id %in% mhc15) | (exon == "brd2" & id %in% brd15))

#get datasets with only classical or classical/intermediate alleles
#in both cases full BRD data is retained
#joins class_nonclass status to gen
gen_c_n <- gen %>% left_join(class_nonclass, by = "label")
#BRD + MHC with all 4 canonical docking aa in ex3 and 1 in ex2
gen_class <- gen_c_n %>% filter(!grepl("nonclass|intermed", class_nonclass)) %>% select(-class_nonclass)
#BRD + MHC with at least 3 of 4 of canonical docking aa in ex3 and 1 in ex2
gen_class_intermed <- gen_c_n %>% filter(!grepl("nonclass", class_nonclass)) %>% select(-class_nonclass)


#Calculate diversity and save reslts ####
#for (dataset in c("gen", "gen15")){
for (dataset in c("gen", "gen15", "gen_class", "gen_class_intermed")) {
  g <- eval(parse(text = dataset))
  suff <- case_when(
    dataset == "gen" ~ "",
    dataset == "gen15" ~ "gen15",
    dataset == "gen_class" ~ "classical",
    dataset == "gen_class_intermed" ~ "classical_and_intermediate"
  )
  a <- lapply(sp_abr1, get_div_MHCBRD, fr = g)
  sp <- lapply(a, "[[", "sp_div")
  ind <- lapply(a, "[[", "ind_div")
  sp_all  <- bind_rows(sp) %>% arrange(species, q, exon, seq_type)
  ind_all  <- bind_rows(ind) %>% arrange(species, id, q, exon, seq_type)
  saveRDS(sp_all, paste0("All_sp_div_MHC_BRD", suff, ".rds"))
  saveRDS(ind_all, paste0("All_sp_indivdual_alphadiv_MHD_BRD", suff, ".rds"))
}


##Df aligned full codon seq of functional alleles ####
#gen_funct <- gen %>% filter(funct == "y") %>% mutate(len = nchar(str_replace_all(seq_codon, "-", ""))) %>%
#  select(genus, species, exon, label, seq_codon, len) %>% distinct() %>% arrange(genus, species, label)
#saveRDS(gen_funct, "Funct_alleles_df_BRD2_MHC_ex2_ex3_all_genera.rds")

#Fasta alignments of functional alleles ####
# #separate file for each species/exon combination
# gen_funct_MHC <- gen %>% filter(funct == "y" & (exon == "2" | exon == "3")) %>% mutate(gensp = paste(genus, species, sep = "_")) %>%
#   select(gensp, exon, label, seq_codon) %>% distinct() %>% arrange(gensp, label)
# saveRDS(gen_funct_MHC, "Funct_alleles_MHC_ex2_ex3_all_genera.rds")
# species <- gen_funct_MHC %>% pull(gensp) %>% unique()
# for(s in species){
#   for(ex in c("2", "3")){
#     sq <- gen_funct_MHC %>% filter(exon == ex & gensp == s) %>% select(label, seq_codon) %>% as.data.frame()
#     dataframe2fas(sq, file = paste0(s, "_ex", ex, ".fas"))
#   }
# }



# #Various checks ####
# #from info
# al_inf <- al %>% filter(grepl('_ex[23]_', label))
# all_nam_inf <- al_inf %>% pull(label)
# f_nam_inf <- al_inf %>% filter(funct == "y") %>% pull(label)
# nf_nam_inf <- al_inf %>% filter(funct == "n") %>% pull(label)
# #functional from alignment
# f_nam_aln <- c(names(ex2_fun_aln), names(ex3_fun_aln))
# #from genotypes
# all_nam_gen <- gen %>% filter(grepl('_ex[23]_', label)) %>% pull(label) %>% unique()
# 
# ##check for alleles without names
# noname <- gen %>% filter(is.na(label))
# 
# ##check allele sets match in different files
# #setdiff(f_nam_aln, f_nam_inf)
# #setdiff(f_nam_inf, f_nam_aln)
# #setdiff(all_nam_inf, all_nam_gen)
# #setdiff(all_nam_gen, all_nam_inf)
# #nonf <- setdiff(all_nam_gen, f_nam_aln)
# #setequal(nonf, nf_nam_inf)


# #Basic stats per species ####
# to_stat <- gen_c_n %>% select(-c(seq, n_reads, subspecies, country)) %>%
#    group_by(id, exon) %>% mutate(ind_n_all = n_distinct(label),
#                                  ind_n_funct_all = n_distinct(label[funct == "y"]),
#                                  ind_n_cons_anch_all = n_distinct(label[class_nonclass == "class"])) %>%
#    ungroup() %>% left_join(fam_gen, by = c("genus", "family"))
# stat1 <- to_stat %>% group_by(family, genus, species, exon, sp_abr1) %>%
#   summarise(N_populations = n_distinct(locality),
#             N_individuals = n_distinct(id),
#             N_alleles = n_distinct(label),
#             N_funct_alleles = n_distinct(label[funct == "y"]),
#             Fr_cons_anch_alleles = n_distinct(label[class_nonclass == "class"])/N_funct_alleles
#             )
# stat2 <- to_stat %>% select(-c(label, funct, locality)) %>%
#   distinct(id, exon, .keep_all = TRUE) %>% group_by(family, genus, species, exon, sp_abr1) %>%
#   summarise(mean_N_alleles_per_ind = mean(ind_n_all),
#             mean_N_funct_alleles_per_ind = mean(ind_n_funct_all),
#             mean_N_cons_anch_alleles_per_ind = mean(ind_n_cons_anch_all),
#             min_N_alleles_per_ind = min(ind_n_all),
#             min_N_funct_alleles_per_ind = min(ind_n_funct_all),
#             min_N_cons_anch_alleles_per_ind = min(ind_n_cons_anch_all),
#             max_N_alleles_per_ind = max(ind_n_all),
#             max_N_funct_alleles_per_ind = max(ind_n_funct_all),
#             max_N_cons_anch_alleles_per_ind = max(ind_n_cons_anch_all),
#             mean_cov = mean(cov))
# 
# Table_S3 <- left_join(stat1, stat2, by = c("family", "genus", "species", "exon", "sp_abr1")) %>% 
#   filter(exon != "brd2") %>% select(-c(sp_abr1, mean_cov)) %>% arrange(exon, family, genus, species)
# 
# write_tsv(Table_S3, "Table_S3.txt")

# #calculate min and max number of allelels for ss 15
# temp <- to_stat %>% select(id, exon, sp_abr1, ind_n_all, ind_n_funct_all) %>%
#   distinct(id, exon, .keep_all = TRUE) %>% group_by(sp_abr1, exon)
# 
# r <- vector("list", 1000)
# for (i in c(1:1000)) {
#   s <- temp %>% sample_n(15) %>% summarise(min_ss15 = min(ind_n_all),
#                                             max_ss15 = max(ind_n_all),
#                                             min_funct_ss15 = min(ind_n_funct_all),
#                                             max_funct_ss15 = max(ind_n_funct_all))
#   r[[i]] <- s
# }
# 
# 
# #generate raw sampling data for PM
# temp1 <- temp %>% filter(exon != "brd2") %>% left_join(taxonomy, by = "sp_abr1") %>% 
#   select(genus, species, sp_abr1, id, exon, ind_n_all, ind_n_funct_all) 
# 
# r_raw <- vector("list", 1000)
# for (i in c(1:1000)) {
#   s_raw <- temp1 %>% sample_n(15) %>% mutate(rep = i)
#   r_raw[[i]] <- s_raw
# }
# 
# r_raw_df <- bind_rows(r_raw) %>% ungroup() %>% select(-c(sp_abr1))
# 
# saveRDS(r_raw_df, "MHC_n_alleles_15ind_resampling.rds")
# 
# stat3 <- bind_rows(r) %>% group_by(sp_abr1, exon) %>% 
#   summarise(min_15ind = mean(min_ss15),
#             max_15ind = mean(max_ss15),
#             min_funct_15ind = mean(min_funct_ss15),
#             max_funct_15ind = mean(max_funct_ss15))
# 
# 
# stat123 <- left_join(stat1, stat2, by = c("family", "genus", "species", "exon", "sp_abr1")) %>% 
#   left_join(stat3, by = c("sp_abr1", "exon"))
# write_lnx_head(stat123, "Basic_statistics_MHC_BRD.txt")


# #generate trees for genera
# genera_trees <- lapply(genera, genus_trees)
# 
# #plot trees for genera to single pdf
# pdf(file = "Trees_ex2_ex3_genera.pdf", width = 10, height = 8, onefile = TRUE)
# for(i in seq_along(genera_trees)){
#   plot(genera_trees[[i]])
# }
# dev.off()
# 
# #check how many alleles are shared between species/genera
# nrow(distinct(gen, label, genus, species))
# nrow(distinct(gen, label, genus))
# gen %>% pull(label) %>% unique() %>% length()

# #constructs and annotates with taxonomy bionj trees of all functional alleles
# #as well as trees of fuctional alleles in Plethodontidae and Salamandridae
# #returns them as list with named elements
# #IMPORTANT - if the trees exist as .rds, this code just loads them into the list,
# 
# d <- gen %>% filter(funct == "y") %>% distinct(label, exon, genus) %>% left_join(fam_gen, by = "genus")
# ts <- vector("list")
# for(ex in c(2, 3)){
#   alg <- get(paste0("ex", ex, "_fun_aln"))
#   ann <- d %>% filter(exon == ex) %>% select(label, family, genus)
#   t_ann_filename <- paste0("All_func_ex", ex, "_annotated_tree.rds")
#   if(file.exists(t_ann_filename)){
#     t_ann <- readRDS(t_ann_filename)
#   } else {
#     JC <- dist.dna(alg, model = "JC", pairwise.deletion = TRUE)
#     tree <- bionj(JC)
#     tree <- midpoint(tree)
#     t_ann <- full_join(tree, ann, by = "label")
#     saveRDS(t_ann, t_ann_filename)
#   }
#   fun_al_fam <- sapply(families, function(x){
#     ann %>% filter(family == x) %>% pull(label) %>% as.character()},
#     simplify = FALSE, USE.NAMES = TRUE)
#   fun_al_gen <- sapply(genera, function(x){
#     ann %>% filter(genus == x) %>% pull(label) %>% as.character()},
#     simplify = FALSE, USE.NAMES = TRUE)
#   t_ann_fam <- groupOTU(t_ann, fun_al_fam, "Family")
#   aaa <- attributes(t_ann_fam@phylo)$Family
#   aaa[aaa == 0] <- "Proteidae"
#   aaa <- droplevels(aaa)
#   attributes(t_ann_fam@phylo)$Family <- aaa
#   pleth <- drop.tip(t_ann, setdiff(t_ann@phylo$tip.label, fun_al_fam[["Plethodontidae"]]))
#   pleth_gen <- groupOTU(pleth, fun_al_gen, "Genus")
#   sal <- drop.tip(t_ann, setdiff(t_ann@phylo$tip.label, fun_al_fam[["Salamandridae"]]))
#   sal_gen <- groupOTU(sal, fun_al_gen, "Genus")
#   ts[[paste0("All_ex",ex)]] <- t_ann_fam
#   ts[[paste0("Plethodontidae_ex",ex)]] <- pleth_gen
#   ts[[paste0("Salamandridae_ex",ex)]] <- sal_gen
# }
# 
# 
# #Plot trees
# #parameters had to be adjusted manually for each tree
# fam_col <- c("Ambystomatidae" = "blueviolet",
#              "Cryptobranchidae" = "chocolate",
#              "Hynobiidae" = "dodgerblue",
#              "Plethodontidae" = "firebrick2",
#              "Proteidae" = "darkgray",
#              "Salamandridae" = "goldenrod")
# 
# t2_fam <- ggtree(ts[["All_ex2"]], layout = "circular", size = 0.5, aes(color = Family)) + xlim(-0.3, NA) +
#   geom_treescale(x = -0.3, width = 0.2, offset = 35, fontsize = 3, linesize = 0.7, color='grey30') +
#   theme(legend.key.size = unit(0.5, "cm"),
#         legend.title = element_blank(),
#         legend.position = c(0.9,0.35),
#         legend.text = element_text(size = 11),
#         legend.box.spacing = unit(-2.8, "cm"),
#         legend.spacing.y = unit(0.2, "cm"),
#         plot.margin = margin(-1, -2, -3.5, -3, "cm")) +
#   scale_color_manual(values = fam_col) +
#   guides(color = guide_legend(override.aes = list(size = 1)))
# t2_fam
# 
# t3_fam <- ggtree(ts[["All_ex3"]], layout = "circular", size = 0.5, aes(color = Family)) + xlim(-0.37, NA) +
#   geom_treescale(x = -0.3, width = 0.2, offset = 35, fontsize = 3, linesize = 0.7, color='grey30') +
#   theme(legend.key.size = unit(0.5, "cm"),
#         legend.title = element_blank(),
#         legend.position = c(0.85, 0.3),
#         legend.text = element_text(size = 11),
#         legend.box.spacing = unit(-2.8, "cm"),
#         legend.spacing.y = unit(0.2, "cm"),
#         plot.margin = margin(0, -2.5, -2, -2.5, "cm")) +
#   scale_color_manual(values = fam_col) +
#   guides(color = guide_legend(override.aes = list(size = 1.5)))
# t3_fam
# 
# tree_ex23 <- ggarrange(t2_fam, t3_fam,
#                        labels = c("exon 2", "exon 3"),
#                        ncol = 1, nrow = 2,
#                        font.label = list(size = 14, color = "black", face = "bold", family = NULL),
#                        hjust = -1,
#                        vjust = 5,
#                        legend = "bottom",
#                        heights = c(1, 1.1),
#                        common.legend = TRUE)
# 
# tree_ex23
# ggsave("Fig_3.pdf", width = 6, height = 10)
# 
# # 
# # t2_plet <- ggtree(ts[["Plethodontidae_ex2"]], layout = "circular", size = 0.4, aes(color = Genus)) + xlim(-0.15, NA) +
# #   geom_treescale(x = -0.15, width = 0.1, offset = 20, fontsize = 3, linesize = 0.7, color='grey30') +
# #   theme(legend.key.size = unit(0.5, "cm"),
# #         legend.position = c(0.8,0.35),
# #         legend.text = element_text(size = 9, face = "italic"),
# #         legend.box.spacing = unit(-2.8, "cm"),
# #         legend.spacing.y = unit(0.2, "cm"),
# #         plot.margin = margin(-3, -3, -5, -6, "cm")) +
# #   guides(color = guide_legend(override.aes = list(size = 1))) 
# # t2_plet
# # 
# # t3_plet <- ggtree(ts[["Plethodontidae_ex3"]], layout = "circular", size = 0.4, aes(color = Genus)) + xlim(-0.15, NA) +
# #   geom_treescale(x = -0.15, width = 0.1, offset = 20, fontsize = 3, linesize = 0.7, color='grey30') +
# #   theme(legend.key.size = unit(0.5, "cm"),
# #         legend.position = c(0.9,0.35),
# #         legend.text = element_text(size = 9, face = "italic"),
# #         legend.box.spacing = unit(-2.8, "cm"),
# #         legend.spacing.y = unit(0.2, "cm"),
# #         plot.margin = margin(-2.5, -1, -2, -3, "cm")) +
# #   guides(color = guide_legend(override.aes = list(size = 1))) 
# # t3_plet
# # 
# # plet_ex23 <- ggarrange(t2_plet, t3_plet,
# #                        labels = c("ex2", "ex3"),
# #                        ncol = 2, nrow = 1,
# #                        font.label = list(size = 11, color = "black", face = "bold", family = NULL),
# #                        hjust = -2,
# #                        vjust = 10,
# #                        legend = "right",
# #                        common.legend = TRUE)
# # 
# # 
# # ggsave("Plethodontids_functional_alleles_ex2_ex3_tree.pdf", width = 15, height = 8)
# # 
# # t2_sal <- ggtree(ts[["Salamandridae_ex2"]], layout = "circular", size = 0.4, aes(color = Genus)) + xlim(-0.2, NA) +
# #   geom_treescale(x = -0.15, width = 0.1, offset = 20, fontsize = 3, linesize = 0.7, color='grey30') +
# #   theme(legend.key.size = unit(0.5, "cm"),
# #         legend.position = c(0.87,0.35),
# #         legend.text = element_text(size = 9, face = "italic"),
# #         legend.box.spacing = unit(-2.8, "cm"),
# #         legend.spacing.y = unit(0.2, "cm"),
# #         plot.margin = margin(-7, 0, -8, -3.5, "cm")) +
# #   guides(color = guide_legend(override.aes = list(size = 1))) 
# # t2_sal
# # 
# # t3_sal <- ggtree(ts[["Salamandridae_ex3"]], layout = "circular", size = 0.4, aes(color = Genus)) + xlim(-0.12, NA) +
# #   geom_treescale(x = -0.15, width = 0.1, offset = 20, fontsize = 3, linesize = 0.7, color='grey30') +
# #   theme(legend.key.size = unit(0.5, "cm"),
# #         legend.position = c(0.9,0.35),
# #         legend.text = element_text(size = 9, face = "italic"),
# #         legend.box.spacing = unit(-2.8, "cm"),
# #         legend.spacing.y = unit(0.2, "cm"),
# #         plot.margin = margin(-2.5, -1, -3.5, -4.5, "cm")) +
# #   guides(color = guide_legend(override.aes = list(size = 1))) 
# # t3_sal
# # 
# # sal_ex23 <- ggarrange(t2_sal, t3_sal,
# #                        labels = c("ex2", "ex3"),
# #                        ncol = 2, nrow = 1,
# #                        font.label = list(size = 11, color = "black", face = "bold", family = NULL),
# #                        hjust = -2,
# #                        vjust = 10,
# #                        legend = "right",
# #                        common.legend = TRUE)
# # 
# # ggsave("Salamandrids_functional_alleles_ex2_ex3_tree.pdf", width = 15, height = 8)
# # 
