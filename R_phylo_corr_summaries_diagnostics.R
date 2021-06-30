#Load libraries and functions####
library(ggpubr)
library(caper)
library(tidyverse)

#load common functions
source("https://raw.githubusercontent.com/wbabik/R_functions/main/R_functions_WB_01.R")

#takes df with PGLS results
#gets table with model characteristics
#parameters, SE and p-val, lambda
#df, adjR2, F, Fdf1, Fdf2, Fp-val
#--------------------
get_mod_characteristics <- function(pgls_file) {
  pg <- readRDS(pgls_file)
  corename <- str_replace(pgls_file, ".rds", "")
  #only full three variable models
  to_diag <- pg %>%  
    select(x, y, cov, seq_type, div_type, q, mod)
  pgg <- pg %>% select(-mod)
  #get table with model characteristics
  #parameters, SE and p-val, lambda
  #df, adjR2, F, Fdf1, Fdf2, Fp-val
  #--------------------
  mchar <- vector("list", nrow(to_diag))
  for (i in 1:nrow(to_diag)){
    r <- to_diag[i, ]
    if (length(r$mod[[1]]) == 2) ms <- summary(r$mod[[1]][["main"]]) else ms <- summary(r$mod[[1]])
    #ms <- summary(r$mod[[1]][["main"]])
    params <- as.data.frame(ms$coefficients)
    params$param <- rownames(params)
    rownames(params) <- NULL
    lambda <- ms$param[["lambda"]]
    R2 <- ms$adj.r.squared
    f <- ms$fstatistic
    fp <- pf(f[[1]], f[[2]], f[[3]], lower.tail = FALSE)
    temp_df <- data.frame(df1 = ms$df[[1]], df2 = ms$df[[2]], lambda = lambda, R2adj = R2,
                          Fstat = f[[1]], Fstat_p_val = fp)
    temp_df <- cbind(params, temp_df)
    i_df <- cbind <- cbind(select(r, -mod), temp_df)
    i_df <- i_df %>% select(x, y, cov, seq_type, div_type, q, 
                            param, Estimate, "Std. Error", "Pr(>|t|)", df1, df2,
                            lambda, R2adj, Fstat, Fstat_p_val)
    colnames(i_df) <- c("predictor", "response", "covariate", "distance", "diversity", "q", 
                        "parameter", "estimate", "SE", "p-val", "df1", "df2", 
                        "lambda", "R2adj", "F", "F_p-val")
    mchar[[i]] <- i_df
  }
  
  model_characteristics <- bind_rows(mchar) %>% 
    arrange(predictor, response, covariate, diversity, distance, q) %>%
    mutate(distance = ifelse(distance == "codon", "dS", distance))
  
  #Technically, for full data this is Table S9
  #and for 15 best
  write_lnx_head(model_characteristics, paste0("model_characteristics_", corename, ".txt"))
}

plot_diag <- function(t, diag, suff) {
  d <- diag %>% filter(seq_type == t[["seq_type"]], div_type == t[["div_type"]], (q==0 | q ==1)) %>% 
    mutate(q = as.factor(q))
    if (t[["div_type"]] == "alpha") d <- d %>% filter(q == 0)
  common_text = paste0(", ", t[["seq_type"]], ", ", t[["div_type"]])
  phyres <- ggplot(d, aes(x = fitted, y = phyres, colour = q)) + geom_point(size = 0.7) +
    geom_smooth() +
    facet_grid(rows = vars(y), scales = "free") +
    theme_bw() + ggtitle(paste0("Phyres vs. fitted", common_text))
  phyresqq <- ggplot(d, aes(sample = phyres, colour = q)) + 
    stat_qq() + stat_qq_line() +
    facet_grid(rows = vars(y), scales = "free") + ggtitle(paste0("Phyres qq", common_text))
  res <-  ggplot(d, aes(x = fitted, y = res, colour = q)) + geom_point(size = 0.7) +
    geom_smooth() +
    facet_grid(rows = vars(y), scales = "free") +
    theme_bw() + ggtitle(paste0("Res vs. fitted", common_text))
  resqq  <- ggplot(d, aes(sample = res, colour = q)) + 
    stat_qq() + stat_qq_line() +
    facet_grid(rows = vars(y), scales = "free") + ggtitle(paste0("Res qq", common_text))
  f <- ggarrange(phyres, res, phyresqq, resqq, ncol = 4)
  ggsave(paste0("Diagnostic_", t[["seq_type"]], "_", t[["div_type"]], suff, ".pdf"), device = "pdf",
         width = 30, height = 20, units = "in")
}

#BODY####

#get text files with model characteristics
pgls_files <- c("PGLS_modeling_results.rds", 
                "PGLS_modeling_results_15_best_covered.rds", 
                "PGLS_modeling_resultsclassical_and_intermediate.rds", 
                "PGLS_modeling_resultsclassical.rds")


lapply(pgls_files, get_mod_characteristics)


pg <- readRDS("PGLS_modeling_results.rds")

to_diag <- pg %>% filter(!is.na(cov)) %>% 
  select(x, y, cov, seq_type, div_type, q, mod)

#df with fitted, obs, res and phyres for all "main" models in to_diag
res <- vector("list", nrow(to_diag))
for (i in 1:nrow(to_diag)) {
  r <- to_diag[i, ]
  m <- r$mod[[1]][["main"]]
  mod_df <- data.frame(fitted = m$fitted, res = m$residuals, phyres = m$phyres, obs = m$y)
  i_df <- cbind(select(r, -mod), mod_df)
  res[[i]] <- i_df
}

diag <- bind_rows(res)

l <- list(c("seq_type" = "AA", "div_type" = "gamma"),
          c("seq_type" = "AA", "div_type" = "alpha"),
          c("seq_type" = "codon", "div_type" = "gamma"),
          c("seq_type" = "codon", "div_type" = "alpha"))

lapply(l, plot_diag, diag, suff = "")


#Prepare Figure_4 APG(MHC-I)####
to_fig <- to_diag %>% filter(y %in% c("APG", "PSMBs", "TAPs", "TAPBP")) %>% 
  filter(seq_type == "AA") %>% filter((div_type == "gamma" & q ==1)| (div_type == "alpha" & q == 0)) %>% 
  distinct(x, y, cov, seq_type, div_type, q, .keep_all = TRUE)

to_fig_plot <- vector("list", nrow(to_fig))
for (i in 1:nrow(to_fig)) {
  r <- to_fig[i, ]
  m <- r$mod[[1]][["main"]]
  print(summary(m))
  pval <- summary(m)$coefficients[2,4]
  print(pval)
  #ax +b
  #if a not significant, then NA
  if (pval > 0.05) {
    a <- NA
    b <- NA
  } else {
    a <- m$model$coef[[2]]
    b <- m$model$coef[[1]]
  }
  print(a)
  print(b)
  d <- as.data.frame(m$x)
  #add NA if not significant
  data <- data.frame(response = m$y, MHC_I = d[[2]], nonAPG = d[[3]], intercept = b, slope = a)
  i_df <- cbind(select(r, -mod), data)
  to_fig_plot[[i]] <- i_df
}

Fig <- bind_rows(to_fig_plot)
FigA <- Fig %>% filter(y == "APG")
FigB <- Fig %>% filter(div_type == "alpha", y %in% c("PSMBs", "TAPs", "TAPBP"))
FigC <- Fig %>% filter(div_type == "gamma", y %in% c("PSMBs", "TAPs", "TAPBP"))
FigB$y  <-  factor(FigB$y, levels = c("PSMBs", "TAPs", "TAPBP"))
FigC$y  <-  factor(FigC$y, levels = c("PSMBs", "TAPs", "TAPBP"))

theme_set(theme_bw(base_size = 16) + theme(plot.title = element_text(hjust = 0.5, size = 20),
                                           axis.title = element_text(size = 16),
                                           panel.grid.major = element_blank(), 
                                           panel.grid.minor = element_blank())) 
update_geom_defaults("point", list(size = 5))
update_geom_defaults("abline", list(size = 2, colour = "grey80"))

fAalpha <- ggplot(filter(FigA, div_type == "alpha"), aes(x = MHC_I, y = response, colour = nonAPG)) + 
  geom_point() + scale_colour_gradient(low= "pink", high = "red4") +
  labs(title = expression(paste("Within-individual (", bold(alpha), italic(", q"), " = 0) diversity"))) +
  xlab("MHC-I") + ylab("APGs average") + theme(axis.title.x = element_text(face = "italic"))

fAalpha
  
fAgamma <- ggplot(filter(FigA, div_type == "gamma"), aes(x = MHC_I, y = response, colour = nonAPG)) + 
  geom_point() +  scale_colour_gradient(low = "#56B1F7", high =  "#132B43") +
  geom_abline(aes(intercept = intercept, slope = slope)) + 
  labs(title = expression(paste("Species-level (", bold(gamma), italic(", q"), " = 1) diversity"))) +
  xlab("MHC-I") + ylab("APGs average") + theme(axis.title.x = element_text(face = "italic"))
  

fAgamma

fB <- ggplot(FigB, aes(x = MHC_I, y = response, colour = nonAPG)) + geom_point() + 
  geom_abline(aes(intercept = intercept, slope=slope)) + 
  scale_colour_gradient(low = "pink", high =  "red4") +
  facet_grid(cols = vars(y), scales = "free") +
  labs(title = expression(paste("Within-individual (", bold(alpha), italic(", q"), " = 0) diversity: three subcategories of APGs"))) +
  xlab("MHC-I") + ylab("APG") + theme(axis.title.x = element_text(face = "italic"),
                                      legend.position = "none",
                                      strip.text = element_text(size = 18, face = "italic"))
fB

fC <- ggplot(FigC, aes(x = MHC_I, y = response, colour = nonAPG)) + geom_point() + 
  geom_abline(aes(intercept = intercept, slope=slope)) + 
  scale_colour_gradient(low = "#56B1F7", high =  "#132B43") +
               facet_grid(cols = vars(y), scales = "free") +
  labs(title = expression(paste("Species-level (", bold(gamma), italic(", q"), " = 1) diversity: three subcategories of APGs"))) +
  xlab("MHC-I") + ylab("APG") + theme(axis.title.x = element_text(face = "italic"),
                                      legend.position = "none",
                                      strip.text = element_text(size = 18, face = "italic"))
fC

gg <- ggarrange(ggarrange(fAalpha, fAgamma, ncol = 2, labels = "A", font.label = list(size = 20)),
                ggarrange(fB, nrow = 1, labels = "B", font.label = list(size = 20)),
                ggarrange(fC, nrow = 1, labels = "C", font.label = list(size = 20)),
                nrow = 3, heights = c(1.1, 1, 1))
                
ggsave("Fig_4.pdf", width = 13, height = 14, units = "in")
  
