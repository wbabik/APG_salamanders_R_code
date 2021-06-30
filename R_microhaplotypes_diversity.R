library(tidyverse)
library(seqinr)
library(ape)
library(hillR)

# The main function is get_div
# - gets rds produced by R_microhaplotypes.R
# - extracts binary genotypes, sequences of alleles, and calculates diversities for each segment
# 
# The body makes filtering and calculated weighted means per gene
# The effect is:
# - dataframe of per segment diversities for segments, genes and species (saved as .rds and .txt)
# - dataframe of per species per gene diversites (saved as .rds and .txt)
#Inputs:
# - *seg_microhap*.rds produced by R_microhaplotypes.R
#IMPORTANT: this version of the script uses only "raw_MIP_seg" - code dealaing with "black" and "white" in in old/
# - gene_ids.txt - translation table for gene names as in reference files into gene symbols and gene class
#when we have data on codon frame, add functional measures

#load common functions
source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

#raw_MIP_seg - segments based on all MIPs
#raw_MIP_seg_min20cov - segments based on MIPs with median coverage >=20 
#categ is a general purpose suffix, which can for example be used to analyse 15 in 
#requires Grantham_Distance_table.csv in the same folder
#this version calculates both species wide alpha and gamma diversity and alpha diversity for each individual
#it also saves aa sequeneces of all segment haplotypes together with their sequences
get_div <- function(t, categ = "raw_MIP_seg_min20cov", fr = frame, stop.rm = TRUE){
  #gets rds produced by R_microhaplotypes.R
  #extracts binary genotypes, sequences of alleles, and calculates diversities
  #returns list of named elements
  #"sp_div"
  #"aa_seq"
  bin_gen <- function(se, df, categ, N = tot_N){
    #for a segment (not collapsed - for collapsed look at the old code) 
    #filters by PAF and coverage and returns a list with named elements:
    # "DNA", "codon", "AA" - each of these is a list of three elements:
    #1) df with binary encoded genotypes (inds in rows, alleles in cols)
    #2) a single row df with basic statistics of the segment
    #3) dataframe with allele seuences
    # add gene, info etc to df2
    d <- df %>% select(-c(MIP, rank, full_haplo, seg_start_MIP, seg_end_MIP)) %>%
      filter((seg_id == se & polymorphic == TRUE & id_loc_cov >=20 & PAF > 0.1)|
               (seg_id == se & polymorphic == FALSE & id_loc_cov >=3)) %>%
      group_by(species, gene, id, id_cov, seg_start_ref, seg_end_ref, seg_id, 
               polymorphic, seg_hap, len, cds_len, codon_hap, aa_hap, lencodon, lenAA) %>%
      summarize(depth = sum(depth)) %>% ungroup() %>% distinct(id, seg_hap, .keep_all = TRUE) %>% 
      as.data.frame()
    sp <- as.character(d[1, "species"])
    gene <- as.character(d[1, "gene"])
    cds_len <- as.numeric(d[1, "cds_len"])
    #if stop.rm = TRUE then haplotypes with stop codons are filtered out early, 
    #so they don't even enter DNA-based calculations
    if (stop.rm == TRUE) d <- d %>% filter(!grepl("\\*", aa_hap))
    dDNA <- d %>% distinct(id, seg_hap, .keep_all = TRUE)
    dcodon <- d %>% distinct(id, codon_hap, .keep_all = TRUE)
    dAA <- d %>% distinct(id, aa_hap, .keep_all = TRUE) 
    if(nrow(dDNA) == 0){
      summary <- data.frame("species" = sp, "gene" = gene, "segment" = se, "len" = NA, "cds_len_bp" = NA, "n_hap" = NA, 
                            "N_typed" = 0, "fr_typed" = 0, "S" = NA, "S_unamb" = NA, "dmax" = NA)
      r <- list("genotypes" = NULL, "summary" = summary, "seq" = NULL)
      res <- list("DNA" <- r, "codon" <- r, "AA" <- r, "AAGhm" <- r)
    } else {
      res <- NULL
      for(s in c("DNA", "codon", "AA")){
        if(s == "DNA"){
          d <- dDNA
          len <- d$len[1]
          h <- "seg_hap"
        } else if(s == "codon"){
          d <- dcodon %>% filter(!is.na(codon_hap))
          if(nrow(d) > 0){
            len <- d$lencodon[1]
            h <- "codon_hap"
          }
        } else if(s == "AA"){
          #this removes also haplotypes with stop codons from AA alignent
          d <- dAA %>% filter(!is.na(aa_hap))
          if(nrow(d) > 0){
            len <- d$lenAA[1]
            h <- "aa_hap"
          }
        }
        if(nrow(d) == 0){
          summary <- data.frame("species" = sp, "gene" = gene, "segment" = se, "len" = NA, "cds_len_bp" = NA, "n_hap" = NA,
                                "N_typed" = 0, "fr_typed" = 0, "S" = NA, "S_unamb" = NA, "dmax" = NA)
          r <- list("genotypes" = NULL, "summary" = summary, "seq" = NULL)
          res[[s]] <- r
        } else {
          a <- d %>% group_by(across(h)) %>% summarize(tot_depth = sum(depth)) %>% 
            arrange(desc(tot_depth), desc(h)) %>% mutate(allele = paste0("all_", sprintf("%02d",row_number())))
          seqdf <- a %>% select(allele, h) %>% as.data.frame()
          if (s == "AA") seq <- df2prot(seqdf) else seq <- df2DNA(seqdf)
          d_a <- left_join(d, a, by = h) %>% select(id, allele) %>% mutate(presence = 1)
          d_a <- d_a %>% pivot_wider(names_from = allele, values_from = presence, values_fill = list(presence = 0)) %>%
            as.data.frame()
          n_hap <- ncol(d_a) - 1
          n_ind <- nrow(d_a)
          summary <- data.frame("species" = sp, "gene" = gene, "segment" = se, "len" = len,  "cds_len_bp" = cds_len, "n_hap" = n_hap, 
                                "N_typed" = n_ind, "fr_typed" = n_ind/N)
          rownames(summary) <- c()
          r <- list("genotypes" = d_a, "summary" = summary, "seq" = seq)
          res[[s]] <- r
          if(s == "AA") res[["AAGhm"]] <- r
        }
      }
    }
    return(res)
  }
  #takes binary genotypes for DNA, codons or AA - 
  #one element of the list produced by bin_gen
  #vectorized function used to identify the beginning of the first codon in the segment
  #handles also cases when the segment ends after the end of the cds
  c_start <- Vectorize(function(g, fs, fe, ssr, sf = sp_frame){
    if(is.na(fs)){
      if(is.na(fe)){
        cs <- NA
      } else {
        gmin <- sf %>% filter(gene_sym == g) %>% pull(seg_start_ref) %>% min()
        cs <- gmin - ssr + 1
      }
    } else {
      cs <- case_when(fs == 1 ~ 1,
                      fs == 2 ~ 3,
                      fs == 3 ~ 2)
    }
    return(cs)
  })
  #vectorized function used to identify the end of the last codon in the segment
  #handles also cases when the segment ends after the end of the cds
  c_end <- Vectorize(function(g, fs, fe, ser, sf = sp_frame){
    if(is.na(fe)){
      if(is.na(fs)){
        ce <- NA
      } else {
        gmax <- sf %>% filter(gene_sym == g) %>% pull(seg_end_ref) %>% max()
        ce <- -(ser - gmax + 1)
      }
    } else {
      ce <- case_when(fe == 1 ~ -2,
                      fe == 2 ~ -3,
                      fe == 3 ~ -1)
    }
    return(ce)
  })
  #vectorized function calculating the length of the segment that covers cds (not necessarily full codons)
  cds_l <- Vectorize(function(g, fs, fe, ssr, ser, sf = sp_frame){
    if (!is.na(fs) & !is.na(fe)) l = ser - ssr else {
      gmin <- sf %>% filter(gene_sym == g) %>% pull(seg_start_ref) %>% min()
      gmax <- sf %>% filter(gene_sym == g) %>% pull(seg_end_ref) %>% max()
      if (is.na(fs)) {
        if (is.na(fe)) l = 0 else {
          l = ser - gmin
        }
      } else {
        l = gmax - ssr
      }
    }
    return(l)
  })
  #read df created by R_microhaplotypes,
  #replaces Xs with Ns
  print(paste0("Processing ", t))
  gene_ids <- read.table("gene_ids.txt", sep = "\t", header = T, stringsAsFactors = F)
  mh <- readRDS(paste0("rds/", t, "_seg_microhap_", categ, ".rds")) %>% 
    left_join(gene_ids, by = "gene") %>% 
    mutate(seg_hap = str_replace_all(seg_hap, "X", "N")) %>% 
    as.data.frame() 
  tot_N <- mh %>% pull(id) %>% unique() %>% length()
  #gets Grantham_table
  Ghm <- read.csv("Grantham_Distance_table.csv", header = T, row.names = 1,encoding = "UTF-8")
  #gets frame positions for the species
  sp_frame <- fr %>% filter(species == t)
  #gets df with distinct haplotypes and other info
  haps <- mh %>% select(species, gene_sym, seg_id, seg_start_ref, seg_end_ref, seg_hap) %>% distinct()
  #gets df with distinct segments and other info
  seg_starts_ends <- haps %>% select(-seg_hap) %>% distinct()
  #adds frame information to segments
  seg_frame <- seg_starts_ends %>% 
    left_join(select(sp_frame,-c(seg_end_ref, frame_end)), by = c("species", "gene_sym", "seg_start_ref")) %>% 
    left_join(select(sp_frame, -c(seg_start_ref, frame_start)), by = c("species", "gene_sym", "seg_end_ref")) %>%
    select(species, gene_sym, seg_id, frame_start, frame_end)
  #print(seg_frame)
  #joins frame infomation to haplotypes
  #extracts full codons into new variable
  #and their aa sequence into another variable
  codon_haps <- haps %>% left_join(seg_frame, by = c("species", "gene_sym", "seg_id")) %>%
    mutate(codon_start = c_start(gene_sym, frame_start, frame_end, seg_start_ref),
           codon_end = c_end(gene_sym, frame_start, frame_end, seg_end_ref),
           trim_hap = str_sub(seg_hap, codon_start, codon_end),
           trim_len = str_length(trim_hap),
           cds_len = cds_l(gene_sym, frame_start, frame_end, seg_start_ref, seg_end_ref),
           codon_hap = ifelse(trim_len >=3, trim_hap, NA),
           aa_hap = v_translate(codon_hap),
           lencodon = nchar(codon_hap)/3,
           lenAA = nchar(aa_hap))
  #joins codon and aa haplotypes to the rest
  #saveRDS(codon_haps, paste0("codon_haps", t, ".rds"))
  mh <- mh %>% left_join(select(codon_haps, species, gene_sym, seg_id, seg_hap, lencodon, cds_len, codon_hap, lenAA, aa_hap), 
                         by = c("species", "gene_sym", "seg_id", "seg_hap"))
  seg_info <- mh %>% distinct(species, gene_sym, seg_id, seg_start_ref, seg_end_ref, cds_len)
  segs <- unique(mh$seg_id)
  segs <- segs[!is.na(segs)]
  dd <- vector("list", length(segs))
  aa_sequences <- vector("list", length(segs))
  for(i in seq_along(segs)){
    seg <- segs[i]
    print(paste0("segment ", seg))
    bg <- bin_gen(seg, mh, categ)
    if(!is.null(bg[["DNA"]][[1]])){
      #save protein haplotypes with additional info
      if(!is.null(bg[["AA"]][["seq"]])){
        fr_typed <- bg[["AA"]][["summary"]][1, "fr_typed"]
        N_typed <- bg[["AA"]][["summary"]][1, "N_typed"]
        all_counts <- bg[["AA"]][["genotypes"]] %>% 
          pivot_longer(-id, names_to = "label", values_to = "pres_abs" ) %>% 
          group_by(label) %>% summarise(count = sum(pres_abs))
        aa_seq <- bg[["AA"]][["seq"]] %>% bin2df() %>% mutate(species = t, 
                                                 seg_id = seg) %>% select(3, 4, 1, 2) %>% 
          left_join(seg_info, by = c("species", "seg_id")) %>% 
          left_join(all_counts, by = "label") %>% 
          mutate(fr_typed = fr_typed,
                 N_typed = N_typed) %>% 
          select(species, gene_sym, seg_id, seg_start_ref, seg_end_ref, label, seq, count, fr_typed, N_typed)
        aa_sequences[[i]] <- aa_seq
      }
      for(z in c("DNA", "codon", "AA", "AAGhm")){
        if(!is.null(bg[[z]][["genotypes"]])){
          #d will contain two lists: alpha_gamma and individial_alpha
          d <- calc_div(bg[[z]], sp = t, typ = z, Ghm_table = Ghm)
          dd[[i]][[z]] <- d
        }
      }
    }
  }
  #for alpha_gamma
  alpha_gamma <- NULL
  #for individual_alpha
  individual_alpha <- NULL
  for(z in c("DNA", "codon", "AA", "AAGhm")){
    #b will contain two lists: alpha_gamma and individial_alpha
    b <- lapply(dd, "[[", z)
    seq_type <- data.frame(seq_type = z)
    a_g <- lapply(b, "[[", "alpha_gamma") %>% bind_rows()
    a_g <- bind_cols(seq_type, a_g)
    print(a_g)
    i_a <- lapply(b, "[[", "individual_alpha") %>% bind_rows()
    i_a <- bind_cols(seq_type, i_a)
    print(i_a)
    alpha_gamma[[z]] <- a_g
    individual_alpha[[z]] <- i_a
  }
  sp_div <- bind_rows(alpha_gamma)
  ind_div <- bind_rows(individual_alpha)
  saveRDS(sp_div, paste0("out/", t, "_sp_div_", categ, ".rds"))
  saveRDS(ind_div, paste0("out/", t, "_individual_alphadiv_", categ, ".rds"))
  saveRDS(aa_sequences,  paste0("out/", t, "_aa_seq_segments_", categ, ".rds"))
  r <- list("sp_div" = sp_div,
            "ind_div" = ind_div,
            "aa_seq" = aa_sequences)
  return(r)
}

#vectorized function used to translate nt sequence trimmed to full codons
v_translate <- Vectorize(function(x){
  l <- str_split(str_to_lower(x), pattern = "")[[1]]
  if(is.na(l[[1]])){
    return(NA)
  } else {
    t <- paste(translate(l, ambiguous = TRUE), collapse='')
    return(t)}    
})


# BODY # 

#Read frame position of each base in the refseq ####
#this is done only once
#IMPORTANT - coordinates here are 1-based, 
#while segments have 0-based, half open
#mutate deals with that
#plus two frame columns are created for convenience
#pull(frame, gene_sym) %>% unique()

frame <- read_delim("frame/frame_table.txt", delim = " ") %>% filter(!is.na(frame)) %>% 
  mutate(seg_start_ref = Npos - 1,
         seg_end_ref = Npos,
         frame_end = frame) %>% rename(frame_start = frame) %>% 
  select(species, gene_sym, seg_start_ref, seg_end_ref, frame_start, frame_end)

#check whether each frame ends at 3rd codon position
frame_check <- frame %>% group_by(species, gene_sym) %>%
  mutate(max = max(seg_end_ref)) %>%
  filter(seg_end_ref == max)

#Calculate and save cds length for gene/species ####

#calculate cds length for each gene in each species
cdses <- frame %>% select(-c(frame_start, frame_end)) %>% group_by(species, gene_sym) %>%
  summarise(cds_length = n()) %>% as.data.frame()

saveRDS(cdses, "CDS_lengths.rds")
write_lnx_head(cdses, "CDS_lengths.txt")

#Calculate diversities ####


taxa <- c("Amb_tex", "Amb_tig", "And_dav", "Bat_att", "Bat_nig", "Des_fus", "Eur_bis", "Hyd_ita", "Hyd_stri",
          "Hyn_lee", "Hyn_ret", "Hyn_tok", "Ich_alp", "Kar_kor", "Lis_bos", "Lis_hel", "Lis_ita", "Omm_nes", "Omm_oph",
          "Plet_cin", "Pleu_wal", "Prot_ang", "Sal_sal", "Tri_cri", "Tri_dob", "Tri_iva", "Tri_kar", "Tri_mac", "Tri_mar", "Tri_pyg")
 
ds <- c("raw_MIP_seg_min20cov", "raw_MIP_seg_min20cov_15_best_covered")
for (i in ds) {
  x <- lapply(taxa, get_div, stop.rm = TRUE)
  sp_div <- lapply(x, "[[", "sp_div")
  ind_div <- lapply(x, "[[", "ind_div")
  suff <- str_replace(i, "raw_MIP_seg_", "")
  saveRDS(bind_rows(sp_div), paste0("All_sp_div_excluding_stop_", suff, ".rds"))
  saveRDS(bind_rows(ind_div), paste0("All_individual_diversities_excluding_stop_", suff, ".rds"))
  if (suff == "min20cov") {
    y <- lapply(taxa, get_div, stop.rm = FALSE)
    aa_seq <- bind_rows(lapply(y, "[[", "aa_seq"))
    aa_seq_STOP <- aa_seq %>% filter(grepl("\\*", seq))
    sp_div_STOP<- bind_rows(lapply(y, "[[", "sp_div"))
    saveRDS(aa_seq, "All_aa_seg_including_stop_min20cov.rds")
    saveRDS(aa_seq_STOP, "All_aa_seg_WITH_stop_codons_min20cov.rds")
    saveRDS(sp_div_STOP, "All_sp_div_including_stop_min20cov.rds")
  }
}
