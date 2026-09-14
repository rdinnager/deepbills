###############################################################################
## 06_flatbm_geodesic.R — N4b: the decisive guard for the geodesic claim.
## Simulate flat-space BM at ALL nodes, build straight parent->child paths,
## measure geodesic L_G along them. If Phi_path = L_G/L_str STILL grows with
## divergence time on flat-BM data, the "geodesic restoration" is a built-in
## metric artifact. If it stays flat, the real +0.39 reflects genuine structure.
###############################################################################
## Two-repo split 2026-09-14: run from the inner repo root (deepbills/deepbills/, which holds _targets.R).
stopifnot(file.exists("_targets.R"))
suppressMessages({library(tidyverse); library(ape); library(phyf); library(Matrix)})
set.seed(1)

## ---- metric setup (identical to 01) ----
lambda <- 1e-2
rho    <- as.numeric(readRDS("data/v3_bayesian/estimated_rho_16dim_rho_schedule_v3.rds"))
bill_z <- read_csv("data/bills_vae_latent_codes_16dim.csv", show_col_types = FALSE)
cent_df <- bill_z %>% dplyr::select(Species, latent_dim, latent_mean) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_mean)
var_df  <- bill_z %>% dplyr::select(Species, latent_dim, latent_var) %>%
  pivot_wider(names_from = latent_dim, values_from = latent_var)
cent_order <- setdiff(names(cent_df), "Species")
C <- as.matrix(cent_df[, cent_order]); V <- as.matrix(var_df[, cent_order]); invV <- 1/V
bpf <- readRDS("data/bill_pf_16dim.rds")
bill_tips <- bpf %>% dplyr::filter(is_tip) %>% dplyr::select(starts_with("latent_")) %>% as.matrix()
rownames(bill_tips) <- (bpf %>% dplyr::filter(is_tip))$label
latvars <- apply(bill_tips, 2, sd)[cent_order]
metric_G <- function(Z){ m<-nrow(Z); dm<-matrix(0,m,nrow(C))
  for(j in seq_len(ncol(Z))){ dm<-dm+ outer(Z[,j],C[,j],"-")^2 / (outer(Z[,j]^2,V[,j],"+")+latvars[j]) }
  1/(lambda + exp(-dm/rho^2) %*% invV) }
path_LG_Lstr <- function(z0,z1){ s<-seq(0,1,length.out=50)
  Z<-outer(1-s,z0)+outer(s,z1); dz<-Z[2:50,]-Z[1:49,]; zmid<-(Z[2:50,]+Z[1:49,])/2
  c(L_str=sqrt(sum((z1-z0)^2)), L_G=sum(sqrt(rowSums(metric_G(zmid)*dz^2)))) }

## ---- ancestry matrices + branch lengths (from training dataloader) ----
dl <- readRDS("data/v3_bayesian/dataloader_dat_v3.rds")
A_end   <- dl$bill_tree_mat            # root->node (end) membership, rows=edges
A_start <- dl$start_mat                # root->parent membership
blens   <- dl$blens                    # per-node branch length, aligned to rows
edge_ids <- rownames(A_end)
rate_j <- sapply(1:16, function(j){ tree<-pf_as_phylo(bpf)
  Y<-bill_tips[tree$tip.label,]; mean(ape::pic(Y[,j], tree, scaled=TRUE)^2) })

## ---- simulate flat BM increments per node, propagate to node states ----
ncol_nodes <- ncol(A_end)
incr <- sapply(1:16, function(j) rnorm(ncol_nodes, 0, sqrt(rate_j[j] * pmax(blens[colnames(A_end)],1e-8))))
end_states   <- as.matrix(A_end   %*% incr)      # 4040 x 16 (child states)
start_states <- as.matrix(A_start %*% incr)      # 4040 x 16 (parent states)
colnames(end_states)<-cent_order; colnames(start_states)<-cent_order

## ---- geodesic + straight length per synthetic edge ----
cat("computing L_G on", nrow(end_states), "synthetic flat-BM edges ...\n")
res <- matrix(NA, nrow(end_states), 2); t0<-Sys.time()
for(i in seq_len(nrow(end_states))){
  res[i,] <- path_LG_Lstr(start_states[i,], end_states[i,])
  if(i%%1000==0) cat("  ",i,"\n")
}
cat("done", round(difftime(Sys.time(),t0,units="mins"),2),"min\n")

## ---- pairwise Phi vs divergence time on synthetic data ----
tree <- pf_as_phylo(bpf); Ntip<-length(tree$tip.label)
child<-tree$edge[,2]; tip_e<-child<=Ntip
clab<-character(length(child)); clab[tip_e]<-tree$tip.label[child[tip_e]]
clab[!tip_e]<-tree$node.label[child[!tip_e]-Ntip]
Lg<-setNames(res[,2],edge_ids)[clab]; Ls<-setNames(res[,1],edge_ids)[clab]
tG<-tree; tG$edge.length<-ifelse(is.na(Lg),0,Lg)
tS<-tree; tS$edge.length<-ifelse(is.na(Ls),0,Ls)
Dt<-cophenetic(tree); DG<-cophenetic(tG); DS<-cophenetic(tS)
ut<-which(upper.tri(Dt),arr.ind=TRUE); tv<-Dt[ut]
dec<-cut(tv,quantile(tv,seq(0,1,.1)),include.lowest=TRUE)
idx<-unlist(lapply(split(seq_len(nrow(ut)),dec),function(g) sample(g,min(3000,length(g)))))
ut<-ut[idx,]; i<-ut[,1]; j<-ut[,2]; tvv<-Dt[cbind(i,j)]
phi<-DG[cbind(i,j)]/DS[cbind(i,j)]
sl<-function(y,x) coef(lm(log10(y)~log10(x)))[2]
cat(sprintf("\nFLAT-BM synthetic:  change-slope Euclidean-path %.3f | geodesic %.3f\n",
            sl(DS[cbind(i,j)],tvv), sl(DG[cbind(i,j)],tvv)))
cat(sprintf("FLAT-BM synthetic:  Phi_path ~ divergence-time slope = %.3f\n", sl(phi,tvv)))
cat("(REAL data gave Phi_path~time slope +0.39, geodesic change-slope 0.66, euc-path 0.27)\n")
cat("=> if flat-BM Phi slope ~0, the real +0.39 is genuine structure; if ~+0.39, it is a metric artifact.\n")
write_csv(tibble(source="flat_bm_synthetic",
                 change_slope_eucpath=sl(DS[cbind(i,j)],tvv),
                 change_slope_geo=sl(DG[cbind(i,j)],tvv),
                 phi_time_slope=sl(phi,tvv)),
          "output/stasis/flatbm_geodesic_guard.csv")
