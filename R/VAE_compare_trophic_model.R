library(tidyverse)
library(castor)
library(targets)
library(phyf)

avonet <- tar_read(bird_beak_avonet)

trophic_dat <- as.factor(avonet$Trophic.Niche)
trophic_levs <- levels(trophic_dat)
trophic_tips <- as.numeric(trophic_dat) 

names(trophic_tips) <- avonet$label
trophic_tips <- trophic_tips[tree$tip.label]

tree <- pf_as_phylo(avonet)

res <- castor::asr_mk_model(tree, trophic_tips,
                            rate_model = "ARD",
                            reroot = FALSE,
                            root_prior = "flat")

res2 <- castor::asr_mk_model(tree, trophic_tips,
                            rate_model = "ARD",
                            reroot = FALSE,
                            root_prior = "empirical")