###############################################################################
## 03_nodeheight_and_figs.R  — fix E1 node-height test + make figures F1/F2/F3
###############################################################################
setwd("G:/Shared drives/COBL Data/Projects/deepbills")
suppressMessages({library(tidyverse); library(ape); library(phyf); library(MASS)})
select <- dplyr::select
theme_set(theme_minimal(base_size = 12))

el <- read_csv("output/stasis/edge_lengths.csv", show_col_types = FALSE) %>%
  filter(T_e > 1e-3, L_str > 1e-6, is.finite(Phi), Phi > 0)

## ---- E1 node-height test (fixed node/height matching) ----
bpf  <- readRDS("data/bill_pf_16dim.rds")
tree <- pf_as_phylo(bpf)
Ntip <- length(tree$tip.label); Nnode <- tree$Nnode
tips <- bpf %>% filter(is_tip)
Y <- as.matrix(tips %>% select(starts_with("latent_"))); rownames(Y) <- tips$label
Y <- Y[tree$tip.label, ]
pics <- sapply(seq_len(ncol(Y)), function(j) ape::pic(Y[, j], tree, scaled = TRUE))
contrast_mag <- sqrt(rowSums(pics^2))                       # length Nnode, node order (Ntip+1..)
heights_all  <- ape::node.depth.edgelength(tree)
node_h <- heights_all[(Ntip + 1):(Ntip + Nnode)]           # internal-node heights, same order
stopifnot(length(node_h) == length(contrast_mag))
nh_lm  <- lm(contrast_mag ~ node_h)
nh_rlm <- MASS::rlm(contrast_mag ~ node_h, maxit = 200)
cat(sprintf("E1 node-height: lm slope=%.4g (p=%.3g); rlm slope=%.4g\n",
            coef(nh_lm)[2], summary(nh_lm)$coef[2,4], coef(nh_rlm)[2]))
cat("  (negative slope = larger contrasts at deep/old nodes = apparent early burst/slowdown)\n")
write_csv(tibble(node_h, contrast_mag), "output/stasis/nodeheight_points.csv")
write_csv(tibble(estimator="E1_nodeheight", slope_lm=coef(nh_lm)[2],
                 p_lm=summary(nh_lm)$coef[2,4], slope_rlm=coef(nh_rlm)[2]),
          "output/stasis/nodeheight_fit.csv")

## ---- F1: change-vs-time (De Lisle robust form), Euclidean vs geodesic (R3) ----
d3 <- el %>% filter(model == "R3")
f1dat <- bind_rows(
  d3 %>% transmute(T_e, L = L_str, ruler = "Euclidean (straight-line)"),
  d3 %>% transmute(T_e, L = L_G,   ruler = "Riemannian geodesic (L_G)"))
fits1 <- f1dat %>% group_by(ruler) %>%
  summarise(slope = coef(lm(log10(L) ~ log10(T_e)))[2], .groups="drop")
f1 <- ggplot(f1dat, aes(log10(T_e), log10(L), colour = ruler)) +
  geom_point(alpha = 0.12, size = 0.5) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 1) +
  geom_abline(slope = 0.5, intercept = median(log10(f1dat$L)) - 0.5*median(log10(f1dat$T_e)),
              linetype = "dashed", colour = "grey40") +
  annotate("text", x = Inf, y = -Inf, hjust=1.05, vjust=-0.6, size=3, colour="grey40",
           label = "dashed = BM reference (slope 0.5)") +
  labs(title = "F1  Change vs time (per edge, R3 reconstruction)",
       subtitle = sprintf("Euclidean slope %.2f  vs  geodesic slope %.2f  (BM=0.50)",
                          fits1$slope[fits1$ruler=="Euclidean (straight-line)"],
                          fits1$slope[fits1$ruler=="Riemannian geodesic (L_G)"]),
       x = "log10 branch length (Myr)", y = "log10 amount of change", colour = NULL) +
  theme(legend.position = "top")
ggsave("output/stasis/F1_change_vs_time.png", f1, width = 8, height = 5.5, dpi = 140)

## ---- F2: keystone — inflation Phi vs interval (R3) ----
key_te  <- coef(lm(log(Phi) ~ log(T_e),      data = d3))[2]
key_age <- coef(lm(log(Phi) ~ log(node_age), data = d3))[2]
f2 <- ggplot(d3, aes(log(T_e), log(Phi))) +
  geom_point(alpha = 0.12, size = 0.5, colour = "#3b6") +
  geom_smooth(method = "lm", se = TRUE, colour = "#164", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  labs(title = "F2  Keystone: does Euclidean underestimation grow with interval?",
       subtitle = sprintf("slope log(Phi)~log(T_e) = %.3f  (Russell's hypothesis predicts > 0; = beta_geo - beta_euc)", key_te),
       x = "log branch length (Myr)", y = "log inflation  Phi = L_G / L_str")
ggsave("output/stasis/F2_keystone_phi.png", f2, width = 8, height = 5.5, dpi = 140)

## ---- F3: rate-vs-interval per model/ruler (E3) ----
f3dat <- el %>% filter(model %in% c("R1","R3")) %>%
  transmute(model, T_e,
            `Euclidean` = rate_str, `Geodesic` = rate_G) %>%
  pivot_longer(c(Euclidean, Geodesic), names_to="ruler", values_to="rate") %>%
  filter(model == "R3" | ruler == "Euclidean")   # R1 only has meaningful Euclidean
f3 <- ggplot(f3dat, aes(log10(T_e), log10(rate), colour = interaction(model, ruler))) +
  geom_smooth(method="lm", se=FALSE, linewidth=1) +
  labs(title="F3  Rate vs interval (blunderbuss slope beta)",
       subtitle="more-negative = stronger apparent slowdown",
       x="log10 branch length (Myr)", y="log10 rate", colour=NULL) +
  theme(legend.position="top")
ggsave("output/stasis/F3_rate_vs_interval.png", f3, width = 8, height = 5.5, dpi = 140)

## ---- E1 figure ----
f_nh <- ggplot(tibble(node_h, contrast_mag), aes(node_h, contrast_mag)) +
  geom_point(alpha = 0.25, size = 0.6) +
  geom_smooth(method = "lm", se = TRUE, colour = "firebrick") +
  labs(title = "E1  Node-height test (observed tip latents, model-free)",
       subtitle = sprintf("slope = %.3g (neg = apparent early burst); root=0, tips=%.0f Myr",
                          coef(nh_lm)[2], max(heights_all)),
       x = "node height above root (Myr)", y = "|standardized contrast| (16-dim)")
ggsave("output/stasis/E1_nodeheight.png", f_nh, width = 8, height = 5.5, dpi = 140)

cat("\nFigures written to output/stasis/: F1, F2, F3, E1\n")
cat(sprintf("keystone slope T_e=%.3f  node_age=%.3f\n", key_te, key_age))
