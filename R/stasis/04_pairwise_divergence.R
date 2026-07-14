###############################################################################
## 04_pairwise_divergence.R — Uyeda-style pairwise divergence-vs-time (E3b/E4)
##
## The proper multi-timescale test across deep splits (per-branch misses this).
##   time_ij      = patristic TIME between tips i,j (real tree; = 2*(H - h_MRCA))
##   euc_ij       = ||z_i - z_j||  (OBSERVED tip latents; MODEL-FREE)
##   euc_path_ij  = summed L_str along path (Euclidean patristic)
##   geo_ij       = summed L_G   along path (geodesic patristic; model-dependent)
## Regress log(rate)~log(time) (blunderbuss beta) and log(change)~log(time)
## (BM ref 0.5), Euclidean vs geodesic ruler.
###############################################################################
setwd("G:/Shared drives/COBL Data/Projects/deepbills")
suppressMessages({library(tidyverse); library(ape); library(phyf)})
select <- dplyr::select
set.seed(1)

bpf  <- readRDS("data/bill_pf_16dim.rds")
tree <- pf_as_phylo(bpf)
Ntip <- length(tree$tip.label)

## observed tip latents (model-free)
tips <- bpf %>% filter(is_tip)
Y <- as.matrix(tips %>% select(starts_with("latent_"))); rownames(Y) <- tips$label
Y <- Y[tree$tip.label, ]

## per-edge lengths (R3) keyed by child-node label
el <- read_csv("output/stasis/edge_lengths.csv", show_col_types = FALSE) %>% filter(model=="R3")
Lg   <- setNames(el$L_G,   el$edge)
Lstr <- setNames(el$L_str, el$edge)

## child label for each tree edge (tips -> tip.label, internal -> node.label)
child <- tree$edge[, 2]
child_lab <- character(length(child))
tip_edge  <- child <= Ntip
child_lab[tip_edge]  <- tree$tip.label[child[tip_edge]]
child_lab[!tip_edge] <- tree$node.label[child[!tip_edge] - Ntip]
te_time <- tree$edge.length
te_Lg   <- Lg[child_lab];   te_Lstr <- Lstr[child_lab]
cat("edges matched: L_G", sum(!is.na(te_Lg)), "/", length(te_Lg), "\n")

## pseudo-trees whose branch lengths are the change measures -> cophenetic = patristic sum
tree_time <- tree
tree_Lg   <- tree; tree_Lg$edge.length   <- ifelse(is.na(te_Lg), 0, te_Lg)
tree_Lstr <- tree; tree_Lstr$edge.length <- ifelse(is.na(te_Lstr),0, te_Lstr)

cat("computing cophenetic matrices (2021 tips) ...\n")
D_time <- cophenetic(tree_time)            # patristic TIME
D_geo  <- cophenetic(tree_Lg)              # summed L_G
D_pstr <- cophenetic(tree_Lstr)            # summed L_str along path

## sample tip pairs across the full time range (stratify by time decile for coverage)
ut <- which(upper.tri(D_time), arr.ind = TRUE)
tvals <- D_time[ut]
dec <- cut(tvals, quantile(tvals, seq(0,1,.1), na.rm=TRUE), include.lowest=TRUE)
idx <- unlist(lapply(split(seq_len(nrow(ut)), dec), function(g) sample(g, min(3000, length(g)))))
ut <- ut[idx, ];
i <- ut[,1]; j <- ut[,2]
lab_i <- rownames(D_time)[i]; lab_j <- colnames(D_time)[j]

pw <- tibble(
  time   = D_time[cbind(i,j)],
  euc_obs = sqrt(rowSums((Y[lab_i,,drop=FALSE] - Y[lab_j,,drop=FALSE])^2)),  # model-free
  euc_path = D_pstr[cbind(i,j)],
  geo    = D_geo[cbind(i,j)]
) %>% filter(time > 0, euc_obs > 0, geo > 0) %>%
  mutate(rate_euc_obs = euc_obs/time, rate_euc_path = euc_path/time, rate_geo = geo/time,
         Phi_path = geo/euc_path)
cat("n pairs:", nrow(pw), " time range:", round(range(pw$time),1), "Myr\n")
write_csv(pw, "output/stasis/pairwise.csv")

sl <- function(y,x){ m<-lm(log10(y)~log10(x)); c(slope=coef(m)[2], p=summary(m)$coef[2,4]) }
cat("\n=== PAIRWISE blunderbuss beta = slope log10(rate)~log10(time) ===\n")
cat(sprintf("  Euclidean (observed tips, model-free): beta = %.3f\n", sl(pw$rate_euc_obs, pw$time)[1]))
cat(sprintf("  Euclidean (summed straight path)     : beta = %.3f\n", sl(pw$rate_euc_path, pw$time)[1]))
cat(sprintf("  Geodesic  (summed L_G path)          : beta = %.3f\n", sl(pw$rate_geo, pw$time)[1]))
cat("\n=== change-vs-time slope (BM ref 0.5) ===\n")
cat(sprintf("  Euclidean observed: %.3f | Euclidean path: %.3f | Geodesic: %.3f\n",
            sl(pw$euc_obs,pw$time)[1], sl(pw$euc_path,pw$time)[1], sl(pw$geo,pw$time)[1]))
cat(sprintf("\n  Phi_path ~ time slope (log-log): %.3f  (hypothesis predicts > 0)\n",
            sl(pw$Phi_path, pw$time)[1]))

## figure: pairwise blunderbuss (rate vs time), observed-Euclidean vs geodesic
theme_set(theme_minimal(base_size=12))
fp <- pw %>% transmute(time, Euclidean=rate_euc_obs, Geodesic=rate_geo) %>%
  pivot_longer(c(Euclidean,Geodesic), names_to="ruler", values_to="rate") %>%
  ggplot(aes(log10(time), log10(rate), colour=ruler)) +
  geom_point(alpha=.06, size=.5) + geom_smooth(method="lm", se=FALSE, linewidth=1) +
  labs(title="F4  Pairwise divergence: rate vs divergence time (blunderbuss)",
       subtitle=sprintf("beta  Euclidean(obs) %.2f  vs  geodesic %.2f  (more neg = stronger apparent slowdown)",
                        sl(pw$rate_euc_obs,pw$time)[1], sl(pw$rate_geo,pw$time)[1]),
       x="log10 divergence time (Myr)", y="log10 rate", colour=NULL) +
  theme(legend.position="top")
ggsave("output/stasis/F4_pairwise_blunderbuss.png", fp, width=8, height=5.5, dpi=140)
cat("\nwrote output/stasis/pairwise.csv + F4_pairwise_blunderbuss.png\n")
