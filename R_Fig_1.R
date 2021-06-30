library(tidyverse)
library(ape)
library(tidytree)
library(ggtree)
tree <- read.nexus("30_salamanders_timetree.tre")
div <- readRDS("Diversities_summary_long_format_15_best_covered.rds")
div <- div %>% filter((gene_sym == "APG" | gene_sym == "MHCI") & q == 1 & seq_type == "AA") %>% 
  select(sp_abr1, gene_sym, PD_alpha, PD_gamma, genus, species, taxon)
div_wide <- div %>% pivot_wider(names_from = gene_sym, values_from = c(PD_alpha, PD_gamma))
#tt <- as_tibble(tree) %>% mutate(var1 = runif(n()),
#                                 var2 = runif(n()))

tt <- as_tibble(tree) %>% left_join(div_wide, by = c("label" = "taxon")) %>% 
  mutate(label = str_replace(label, "_", " "))

ttt <- as.treedata(tt)
ggtree(ttt, size = 1, colour = "grey20") + 
  theme_tree2() + 
  coord_cartesian(clip = 'off') +
  theme(axis.text.x = element_text(face = "bold", size = 11),
        axis.title.x = element_text(hjust = 0.28, size = 13)) +
  geom_tippoint(aes(size = PD_alpha_APG, x = x + 10), pch = 19, colour = "red3") + 
  geom_tippoint(aes(size = PD_alpha_MHCI/5, x = x + 20), pch = 19, colour = "red3") + 
  geom_tippoint(aes(size = PD_gamma_APG, x = x + 30), pch = 19, colour = "navyblue") +
  geom_tippoint(aes(size = PD_gamma_MHCI/5, x = x + 40), pch = 19, colour = "navyblue") +
  geom_tiplab(aes(label = label), offset = 50, fontface = "bold.italic", size = 3.5) +
  scale_x_continuous(breaks=seq(from = 0, to = 200, by = 50), 
                     labels = seq(from = 200, to = 0, by = -50), limits = c(-20, 400)) +
  scale_y_continuous(limits = c(-0.1, 31)) +
  xlab("million years ago") +
  guides(size = FALSE) + theme(plot.margin = unit(c(1.5,0,1,0), "cm")) +
  annotate(geom="text", x=205, y=31, label=expression(paste(alpha, " APGs")),
                                  color="red", angle = 90, hjust =0) +
  annotate(geom="text", x=215, y=31, label=expression(paste(alpha, italic(" MHC-I"))),
           color="red", angle = 90, hjust =0) +
  annotate(geom="text", x=225, y=31, label=expression(paste(gamma, " APGs")),
           color="navyblue", angle = 90, hjust =0) +
  annotate(geom="text", x=235, y=31, label=expression(paste(gamma, italic(" MHC-I"))),
           color="navyblue", angle = 90, hjust =0)
  


ggsave("Fig_1.pdf", width = 8, height = 7, units = "in")
  
