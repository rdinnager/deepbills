## Two-repo split 2026-09-14: run from the inner repo root (deepbills/deepbills/, which holds _targets.R).
stopifnot(file.exists("_targets.R"))
suppressMessages({library(tidyverse); library(phyf); library(ape)})

r3 <- readRDS("data/v3_bayesian/bill_vae_aces_16dim_v3_bayesian.rds")
et <- readRDS("data/bill_edge_trajs_16dim.rds")
r1 <- readRDS("data/bill_w_preds_16dim.rds")

cat("=== edge_trajs start/end/time head ===\n")
print(head(et[,c("start","end","start_time","end_time")]))
cat("start_time range:", range(et$start_time), " end_time range:", range(et$end_time), "\n")
cat("branch length (end-start) summary:\n"); print(summary(et$end_time - et$start_time))

cat("\n=== do et$end match r3$edge (row-aligned)? ===\n")
cat("identical(et$end, r3$edge):", identical(as.character(et$end), as.character(r3$edge)), "\n")
cat("et$end head:", paste(head(et$end),collapse=", "), "\n")
cat("r3$edge head:", paste(head(r3$edge),collapse=", "), "\n")

cat("\n=== r1 time column ===\n")
print(head(data.frame(label=r1$label, is_tip=r1$is_tip, time=r1$time)))
cat("r1 time range:", range(r1$time, na.rm=TRUE), "\n")

cat("\n=== build tree, node depths + branch lengths ===\n")
tree <- pf_as_phylo(readRDS("data/bill_pf_16dim.rds"))
cat("tree tips:", length(tree$tip.label), " nodes:", tree$Nnode, " ultrametric:", is.ultrametric(tree), "\n")
cat("total tree height:", max(node.depth.edgelength(tree)), "\n")
bl <- tree$edge.length
cat("branch length summary (tree):\n"); print(summary(bl))
