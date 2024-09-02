####### Load Libraries ###################
library(targets)
library(torch)
library(tidyverse)
library(dagnn) ## This is in-development package: https://github.com/rdinnager/dagnn
library(zeallot)
library(conflicted)
library(phyf)
library(torchopt)

source("R/run_vae.R")

conflicts_prefer(torch::optim_adamw)

####### Load Data ########################

avonet <- tar_read(bird_beak_avonet)

mod <- run_vae(avonet)

options(torch.serialization_version = 2)
torch_save(mod, "data/bill_vae_uncond_v1.to")
