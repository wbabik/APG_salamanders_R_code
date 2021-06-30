library(seqinr)
library(tidyverse)

source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

# FUNCTIONS ####

#takes gene name and dataframe of sequences (names: gene, species, seq)
#extracts gene, converts into AAbin
#aligns with muscle
#saves alignment to file
#generates dataframe with gene name, species name, 
#alignment aa (as column), position in alignment and position in protein
#returns data frame
alg_tru_pos <- function(g, seq = aa_seq) {
  gene <- seq %>% filter(gene == g) %>% select(-gene) %>% df2prot()
  
  #align and save alignment
  gene_alg <- muscle(gene, exec = "muscle3.8.31_i86win32.exe")
  write.FASTA(gene_alg, paste0(g, "_prot.fas"))
  spec <- dimnames(gene_alg)[[1]]
  pos <- function(s){
    u <- c(as.character(gene_alg[s, ]))
    udf <- data.frame(gene_sym = g,
                      species = s,
                      aa = u) %>% 
      mutate(pos_aln = row_number(),
             aa_pos = cumsum(aa != "-"))
    return(udf)
  }
  res_list <- lapply(spec, pos)
  res_df <- bind_rows(res_list)
}


#takes gene symbol, species, segment start and end on reference
#and table with amino acid positions on reference
#produces ; delimited string with aa start and end
#the function is then Vectorized to work with mutate
res_pos <- function(gs, sp, s, e, ft = aa_pos){
  seg <- ft %>% filter(gene_sym == gs, species == sp, seg_start_ref >= s, seg_end_ref <= e)
  res_start <- seg %>% filter(frame == 1) %>% pull(res) %>% min()
  res_end <- seg %>% filter(frame == 3) %>% pull(res) %>% max()
  res <- paste(res_start, res_end, sep = ";")
  return(res)
}
res_pos_v <- Vectorize(res_pos)

# BODY

#Read frame position of each base in the refseq ####
#this is done only once
#IMPORTANT - coordinates here are 1-based, 
#while segments have 0-based, half open
#mutate deals with that
#plus two frame columns are created for convenience
#gets positions within cds
#extracts TAP1 and TAP2

cds <- read_delim("frame/frame_table.txt", delim = " ") %>% filter(!is.na(frame)) %>% 
  mutate(seg_start_ref = Npos - 1,
         seg_end_ref = Npos,
         frame_end = frame) %>% 
  filter(gene_sym == "TAP1" | gene_sym == "TAP2")

#gets list of species and genes
sp <- cds %>% pull(species) %>% unique()
gene <- cds %>% pull(gene_sym) %>% unique()

#for each gene and species extracts nt sequence as vector and translates into protein
#saves to data frame with gene and species info and returns list of data frames
df_list <- vector("list", length(sp)*length(gene))
count <- 0
for (g in gene) {
  for (s in sp) {
    nt <- cds %>% filter(gene_sym == g, species == s) 
    if (nrow(nt) > 0){
      count <- count + 1
      aa <-  nt %>% pull(noNbase) %>% translate(ambiguous = TRUE) %>% paste(collapse = "")
      df <- data.frame(gene = g, species = s, seq = aa)
      df_list[[count]] <- df
    }
  }
}
aa_seq <- bind_rows(df_list)

#adds for each base its aa position
aa_pos <- cds %>% group_by(gene_sym, species) %>% mutate(res = cumsum(frame == 1))

#reads TAP12 aa segment haplotypes
aa_seg <- readRDS("rds/All_species_aa_segments_TAP12.rds")
#ad <- readRDS("Andrias_aa_segments.rds")

#ads start end amino acid positions for each segment
aa_seg_aa_pos <- aa_seg %>% mutate(aa_start_stop = res_pos_v(gs = gene_sym, sp = species, s = seg_start_ref, e = seg_end_ref)) %>% 
  separate(aa_start_stop, c("aa_start", "aa_end"), convert = TRUE) 

#for a single row df creates df with individual aa in rows
aa_pos_df <- function(r) {
  aa <- c(str_split(r[1, "seq"], "", simplify = TRUE))
  aa_pos <- c(r[1, "aa_start"]:r[1, "aa_end"])
  r_df <- cbind(r, aa, aa_pos)
}

#splits aa_seg_aa_pos into list of rows
aa_seg_aa_pos_list <- split(aa_seg_aa_pos, seq(nrow(aa_seg_aa_pos)))

#creates big dataframe with each amino-acid in each haplotype in separate row
#with aa position on reference
d <- lapply(aa_seg_aa_pos_list, aa_pos_df)
aa_df <- bind_rows(d)

#for each aa position calculates the fraction of individuals that have a particular aa
aa_fr_sp <- aa_df %>% mutate(fr = count/N_typed) %>% 
  group_by(species, gene_sym, seg_id, aa_pos, aa) %>% 
  summarise(sum_freq = sum(fr)) %>% 
  mutate(freq_ind = ifelse(sum_freq <= 1, sum_freq, 1)) %>% 
  select(-sum_freq)

#makes aa alignment and
#for each aa in reference provides position in alignment
temp <- lapply(c("TAP1", "TAP2"), alg_tru_pos)

TAP_positions <- bind_rows(temp) %>% select(-aa)

#adds information about position in alignemnt, filters out aa variants present in less than 10% individuals
#this can be filtered by alignment position to generate supplementary table with TAP polymorphism
aa_pol <- aa_fr_sp %>% left_join(TAP_positions, by = c("species", "gene_sym", "aa_pos")) %>% 
  filter(freq_ind >= 0.1) %>% mutate(aa = str_to_upper(aa))

saveRDS(aa_pol, "TAP12_aa_polymorphism.rds")

TAP1_funct <- c(49,57,134,225,228,229,239,256)
TAP2_funct <- c(230,231,232,234,235,258)

Table_S3 <- aa_pol %>% ungroup() %>% 
  filter((gene_sym == "TAP1" & pos_aln %in% TAP1_funct) | (gene_sym == "TAP2" & pos_aln %in% TAP2_funct)) %>%
  select(-seg_id) %>% 
  arrange(gene_sym, pos_aln, species, desc(freq_ind)) %>% select(-c(aa_pos, freq_ind)) %>% 
  group_by(species, gene_sym, pos_aln) %>% mutate(aa_variants = paste(aa, collapse = "/")) %>% 
  select(-aa) %>% distinct() %>% 
  pivot_wider(names_from = species, values_from = aa_variants)

Table_S3[is.na(Table_S3)] <- "?"


