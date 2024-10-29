library(tidyverse)
library(phyf)
library(fibre)
library(Matrix)
library(MCMCglmm)

bird_tree <- pf_as_phylo(bird_beak_codes)
bird_inv <- inverseA(bird_tree, nodes = "tips")
bird_inv_mat <- bird_inv$Ainv

bird_inv_chol <- chol(bird_inv_mat)