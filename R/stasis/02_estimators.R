###############################################################################
## 02_estimators.R  —  stasis rate-vs-time estimators (CPU)
##
## Keystone: log(Phi) ~ log(dt), Phi = L_G/L_str.  Identity: this slope EQUALS
##   beta_geodesic - beta_euclidean (the difference in rate-vs-interval scaling
##   between the geodesic and Euclidean rulers), because
##   log(rate_G) - log(rate_str) = log(L_G/T) - log(L_str/T) = log(Phi).
## So slope > 0  <=>  geodesic ruler flattens the apparent slowdown.
##
## Also: E3 rate-vs-interval beta per ruler/model; E4 change-vs-time slope
## (BM reference = 0.5); E1 Freckleton-Harvey node-height test on observed tips
## (model-free Euclidean signal); N2 numerator-randomization null on keystone.
###############################################################################
setwd("G:/Shared drives/COBL Data/Projects/deepbills")
suppressMessages({library(tidyverse); library(ape); library(phyf); library(MASS)})
select <- dplyr::select
dir.create("output/stasis", showWarnings = FALSE, recursive = TRUE)

el <- read_csv("output/stasis/edge_lengths.csv", show_col_types = FALSE) %>%
  filter(T_e > 1e-3, L_str > 1e-6, is.finite(Phi), Phi > 0)      # drop degenerate edges

lm_slope <- function(y, x, robust = FALSE) {
  d <- tibble(y, x) %>% filter(is.finite(y), is.finite(x))
  m <- if (robust) MASS::rlm(y ~ x, data = d, maxit = 100) else lm(y ~ x, data = d)
  ci <- tryCatch(confint(m)["x", ], error = function(e) c(NA, NA))
  s  <- summary(m)
  p  <- tryCatch(coef(s)["x", 4], error = function(e) NA_real_)
  c(slope = unname(coef(m)["x"]), lo = ci[1], hi = ci[2], p = p, n = nrow(d))
}

results <- list()
add <- function(label, model, estimator, ruler, v)
  results[[length(results)+1]] <<- tibble(label, model, estimator, ruler,
    slope = v["slope"], lo = v["lo"], hi = v["hi"], p = v["p"], n = v["n"])

###############################################################################
## KEYSTONE — does inflation grow with the measurement interval? (R3)
###############################################################################
for (subset in c("all", "internal")) {
  d <- el %>% filter(model == "R3")
  if (subset == "internal") d <- d %>% filter(!is_tip)
  add(paste0("keystone/", subset), "R3", "log(Phi)~log(T_e)", "geo/euc",
      lm_slope(log(d$Phi), log(d$T_e)))
  add(paste0("keystone_metric/", subset), "R3", "log(phi_metric)~log(T_e)", "geo/euc",
      lm_slope(log(d$phi_metric), log(d$T_e)))
  add(paste0("keystone_age/", subset), "R3", "log(Phi)~log(node_age)", "geo/euc",
      lm_slope(log(d$Phi), log(d$node_age)))
}

###############################################################################
## E3 rate-vs-interval (beta) and E4 change-vs-time (slope, BM ref 0.5)
## per ruler, per model.  rate = L/T_e ; change = L.
###############################################################################
rulers <- list(str = c(L="L_str", rate="rate_str"),
               euc = c(L="L_euc", rate="rate_euc"),
               G   = c(L="L_G",   rate="rate_G"))
for (mdl in c("R1", "R3")) {
  d <- el %>% filter(model == mdl)
  for (rn in names(rulers)) {
    Lcol <- rulers[[rn]]["L"]; Rcol <- rulers[[rn]]["rate"]
    add(paste0("E3_beta/", mdl), mdl, "log10(rate)~log10(T_e)", rn,
        lm_slope(log10(d[[Rcol]]), log10(d$T_e)))
    add(paste0("E4_change/", mdl), mdl, "log10(L)~log10(T_e)", rn,
        lm_slope(log10(d[[Lcol]]), log10(d$T_e)))
  }
}

est <- bind_rows(results)
write_csv(est, "output/stasis/estimator_table.csv")
cat("\n=== ESTIMATOR TABLE ===\n"); print(est, n = 40, width = 130)

## Verify the identity beta_geo - beta_euc == keystone slope (R3, all edges)
d <- el %>% filter(model == "R3")
b_euc <- lm_slope(log10(d$rate_str), log10(d$T_e))["slope"]
b_geo <- lm_slope(log10(d$rate_G),   log10(d$T_e))["slope"]
key   <- lm_slope(log(d$Phi), log(d$T_e))["slope"]
cat(sprintf("\nIdentity check: beta_geo - beta_euc = %.4f ; keystone slope = %.4f (should match)\n",
            b_geo - b_euc, key))

###############################################################################
## N2 — numerator-randomization null for the keystone (shuffle T_e vs Phi)
###############################################################################
set.seed(1)
d <- el %>% filter(model == "R3")
obs <- lm_slope(log(d$Phi), log(d$T_e))["slope"]
B <- 1000
null_sl <- replicate(B, lm_slope(log(d$Phi), log(sample(d$T_e)))["slope"])
p_null <- mean(abs(null_sl) >= abs(obs))
cat(sprintf("\nN2 keystone: observed slope %.4f ; null mean %.4f sd %.4f ; p=%.3f\n",
            obs, mean(null_sl), sd(null_sl), p_null))
write_csv(tibble(observed = obs, null_mean = mean(null_sl), null_sd = sd(null_sl),
                 null_lo = quantile(null_sl,.025), null_hi = quantile(null_sl,.975),
                 p = p_null), "output/stasis/keystone_null.csv")

###############################################################################
## E1 — Freckleton-Harvey node-height test on OBSERVED tip latents (model-free)
##   |standardized PIC| (L2 across 16 dims) vs node height above root.
##   Negative slope = early burst / apparent slowdown under the Euclidean ruler.
###############################################################################
tree <- pf_as_phylo(readRDS("data/bill_pf_16dim.rds"))
tips <- readRDS("data/bill_pf_16dim.rds") %>% filter(is_tip)
Y <- as.matrix(tips %>% select(starts_with("latent_")))
rownames(Y) <- tips$label
Y <- Y[tree$tip.label, ]

## per-dim contrasts, node heights
pics <- sapply(seq_len(ncol(Y)), function(j) ape::pic(Y[, j], tree, scaled = TRUE))
contrast_mag <- sqrt(rowSums(pics^2))                 # multivariate |contrast|
node_ids <- as.integer(names(ape::pic(Y[,1], tree, scaled = TRUE)))
heights  <- ape::node.depth.edgelength(tree)          # root=0 .. tips=height
node_h   <- heights[node_ids]
nh <- lm_slope(contrast_mag, node_h)
nh_r <- lm_slope(contrast_mag, node_h, robust = TRUE)
add("E1_nodeheight/tips", "R1-tips(observed)", "|PIC|~node_height", "euclidean", nh)
add("E1_nodeheight_rlm/tips", "R1-tips(observed)", "|PIC|~node_height(rlm)", "euclidean", nh_r)
cat(sprintf("\nE1 node-height (observed latents): slope=%.4g p=%.3g (rlm slope=%.4g)\n",
            nh["slope"], nh["p"], nh_r["slope"]))

est <- bind_rows(results)
write_csv(est, "output/stasis/estimator_table.csv")
cat("\nwrote output/stasis/estimator_table.csv\n")
