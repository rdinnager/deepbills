###############################################################################
## 05_flatbm_control.R — N4 negative control (Euclidean): does flat BM reproduce
## the observed straight-line divergence SATURATION? Simulate BM via VCV+chol.
###############################################################################
setwd("G:/Shared drives/COBL Data/Projects/deepbills")
suppressMessages({library(tidyverse); library(ape); library(phyf)})
set.seed(1)
bpf <- readRDS("data/bill_pf_16dim.rds"); tree <- pf_as_phylo(bpf)
tips <- bpf %>% dplyr::filter(is_tip)
Y <- as.matrix(tips %>% dplyr::select(starts_with("latent_"))); rownames(Y)<-tips$label
Y <- Y[tree$tip.label, ]
n <- length(tree$tip.label)

## per-dim BM rate from contrasts
rates <- sapply(1:16, function(j){ pc<-ape::pic(Y[,j],tree,scaled=TRUE); mean(pc^2) })

## simulate flat multivariate BM tips: sim = t(chol(V)) %*% Z * sigma
cat("building VCV + chol ...\n")
V <- vcv(tree)                          # n x n
Lc <- t(chol(V))                        # lower
sim <- sapply(1:16, function(j) as.numeric(Lc %*% rnorm(n)) * sqrt(rates[j]))
rownames(sim) <- rownames(V)
sim <- sim[tree$tip.label, ]

## reuse cophenetic time + stratified pair sample
Dt <- cophenetic(tree)
ut <- which(upper.tri(Dt), arr.ind=TRUE); tv <- Dt[ut]
dec <- cut(tv, quantile(tv, seq(0,1,.1)), include.lowest=TRUE)
idx <- unlist(lapply(split(seq_len(nrow(ut)), dec), function(g) sample(g, min(3000,length(g)))))
ut <- ut[idx,]; i<-ut[,1]; j<-ut[,2]; li<-rownames(Dt)[i]; lj<-colnames(Dt)[j]; tvv<-Dt[cbind(i,j)]

euc_real <- sqrt(rowSums((Y[li,]-Y[lj,])^2))
euc_sim  <- sqrt(rowSums((sim[li,]-sim[lj,])^2))
sl <- function(y,x) coef(lm(log10(y)~log10(x)))[2]
cat(sprintf("\nREAL    Euclidean change-slope %.3f  (rate beta %.3f)\n", sl(euc_real,tvv), sl(euc_real/tvv,tvv)))
cat(sprintf("FLAT-BM Euclidean change-slope %.3f  (rate beta %.3f)\n", sl(euc_sim,tvv),  sl(euc_sim/tvv,tvv)))
cat("BM reference change-slope = 0.5 ; rate beta = -0.5\n")
cat("\n=> if REAL slope << 0.5 while FLAT-BM ~0.5, the straight-line saturation is genuine\n")
cat("   apparent stasis (bounded/manifold), not a null artifact.\n")
write_csv(tibble(measure=c("real","flat_bm"),
                 change_slope=c(sl(euc_real,tvv), sl(euc_sim,tvv)),
                 rate_beta=c(sl(euc_real/tvv,tvv), sl(euc_sim/tvv,tvv))),
          "output/stasis/flatbm_control.csv")
