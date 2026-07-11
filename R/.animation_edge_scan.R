## Phase 1: Edge Scanning — Find Interesting Edges
## No GPU needed. Reads pre-computed prediction data.
##
## Computes curvature ratio (curved path length / straight-line distance)
## for every edge in the manifold model. High ratios indicate edges where
## evolution "goes around the long way" through latent space.

library(tidyverse)
library(ape)
library(phyf)

source("R/.animation_utils.R")

# ---- Load data ----

message("Loading manifold predictions...")
z_mani <- read_rds("data/evo_predictions_all_edges_v2.rds")

message("Loading tree data...")
bill_edge_trajs <- read_rds("data/bill_edge_trajs_16dim.rds")
trophic_dat <- read_rds("data/trophic_dat.rds")

niche_pal <- load_niche_pal()

# ---- 1A. Curvature ratio per edge ----

message("Computing curvature scores for ", nrow(z_mani), " edges...")

edge_scores <- tibble(
  edge = z_mani$edge,
  is_tip = z_mani$is_tip,
  dist_manifold = map_dbl(z_mani$z_seqs, calc_dist_along, .progress = TRUE),
  straight_dist = map_dbl(z_mani$z_seqs, calc_straight_dist, .progress = TRUE)
) |>
  mutate(
    curvature_ratio = dist_manifold / straight_dist,
    excess_dist = dist_manifold - straight_dist
  )

# Add time info for rate calculation
edge_scores <- edge_scores |>
  mutate(
    time_span = map_dbl(z_mani$time_seqs, ~ .x[length(.x)] - .x[1]),
    rate_manifold = dist_manifold / time_span
  )

# ---- 1B. Trophic transition edges ----

message("Finding trophic transition edges...")

edge_scores <- edge_scores |>
  mutate(
    start_troph = map_chr(z_mani$prtroph, ~ .x[1]),
    end_troph = map_chr(z_mani$prtroph, ~ .x[length(.x)]),
    has_transition = start_troph != end_troph,
    transition = ifelse(has_transition,
                        paste(start_troph, "->", end_troph),
                        NA_character_)
  )

# Add true trophic niche for tips
edge_scores <- edge_scores |>
  left_join(
    trophic_dat |> select(edge = label, true_trophic = trophic_niche),
    by = "edge"
  )

# ---- 1C. Sort and rank ----

edge_scores <- edge_scores |>
  arrange(desc(curvature_ratio)) |>
  mutate(rank = row_number())

# ---- Save results ----

write_csv(edge_scores, "output/edge_scan/edge_curvature_scores.csv")
message("Saved edge scores to output/edge_scan/edge_curvature_scores.csv")

# ---- Summary statistics ----

message("\n=== CURVATURE SUMMARY ===")
message("Total edges: ", nrow(edge_scores))
message("Median curvature ratio: ", round(median(edge_scores$curvature_ratio, na.rm = TRUE), 3))
message("Mean curvature ratio: ", round(mean(edge_scores$curvature_ratio, na.rm = TRUE), 3))
message("Max curvature ratio: ", round(max(edge_scores$curvature_ratio, na.rm = TRUE), 3))
message("Edges with ratio > 1.5: ", sum(edge_scores$curvature_ratio > 1.5, na.rm = TRUE))
message("Edges with ratio > 2.0: ", sum(edge_scores$curvature_ratio > 2.0, na.rm = TRUE))
message("Edges with trophic transitions: ", sum(edge_scores$has_transition, na.rm = TRUE))
message("Edges with transition AND ratio > 1.5: ",
        sum(edge_scores$has_transition & edge_scores$curvature_ratio > 1.5, na.rm = TRUE))

# ---- Top candidates ----

message("\n=== TOP 20 HIGHEST CURVATURE EDGES ===")
top20 <- edge_scores |>
  head(20) |>
  select(rank, edge, is_tip, curvature_ratio, excess_dist,
         dist_manifold, straight_dist, has_transition, transition, true_trophic)
print(top20, n = 20, width = 120)

message("\n=== TOP TROPHIC TRANSITION EDGES (sorted by curvature) ===")
top_transitions <- edge_scores |>
  filter(has_transition) |>
  head(20) |>
  select(rank, edge, is_tip, curvature_ratio, transition, true_trophic)
print(top_transitions, n = 20, width = 120)

# ---- Top edges by absolute excess distance ----
# These have the most visual impact: large absolute detours
# rather than just high ratios on tiny edges

message("\n=== TOP 10 EDGES BY EXCESS DISTANCE ===")
top_excess <- edge_scores |>
  arrange(desc(excess_dist)) |>
  head(10) |>
  select(edge, is_tip, curvature_ratio, excess_dist,
         dist_manifold, straight_dist,
         has_transition, transition)
print(top_excess, n = 10, width = 120)

message("\n=== TOP 10 TIP EDGES BY EXCESS DISTANCE ===")
top_tip_excess <- edge_scores |>
  filter(is_tip) |>
  arrange(desc(excess_dist)) |>
  head(10) |>
  select(edge, curvature_ratio, excess_dist,
         dist_manifold, has_transition, transition)
print(top_tip_excess, n = 10, width = 120)

# Excess distance plot
p_excess <- ggplot(
  edge_scores,
  aes(x = straight_dist, y = excess_dist,
      color = has_transition)
) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed",
             color = "grey50") +
  scale_color_manual(
    values = c("FALSE" = "grey60", "TRUE" = "red"),
    labels = c("No transition", "Trophic transition"),
    name = ""
  ) +
  labs(
    x = "Straight-line Distance (16-dim Euclidean)",
    y = "Excess Distance (manifold - straight)",
    title = "Absolute Excess Path Length for All Edges"
  ) +
  theme_minimal()

ggsave("output/edge_scan/excess_distance.pdf",
       p_excess, width = 8, height = 6)
ggsave("output/edge_scan/excess_distance.png",
       p_excess, width = 8, height = 6, dpi = 150)

# ---- Plots ----

# Histogram of curvature ratios
p_hist <- ggplot(edge_scores, aes(x = curvature_ratio)) +
  geom_histogram(bins = 100, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "red") +
  scale_x_log10() +
  labs(x = "Curvature Ratio (log scale)",
       y = "Count",
       title = "Distribution of Edge Curvature Ratios",
       subtitle = "Ratio = curved path length / straight-line distance; >1 means curved") +
  theme_minimal()

ggsave("output/edge_scan/curvature_histogram.pdf", p_hist, width = 8, height = 5)
ggsave("output/edge_scan/curvature_histogram.png", p_hist, width = 8, height = 5, dpi = 150)

# Curvature ratio vs distance
p_scatter <- ggplot(edge_scores, aes(x = straight_dist, y = curvature_ratio,
                                      color = has_transition)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_y_log10() +
  scale_color_manual(values = c("FALSE" = "grey60", "TRUE" = "red"),
                     labels = c("No transition", "Trophic transition"),
                     name = "") +
  labs(x = "Straight-line Distance (16-dim Euclidean)",
       y = "Curvature Ratio (log scale)",
       title = "Curvature vs Distance for All Edges") +
  theme_minimal()

ggsave("output/edge_scan/curvature_vs_distance.pdf", p_scatter, width = 8, height = 6)
ggsave("output/edge_scan/curvature_vs_distance.png", p_scatter, width = 8, height = 6, dpi = 150)

# ---- Clade-level curvature scores ----

message("\nComputing clade-level curvature scores...")

# Build phylo tree for clade operations
phylo_tree <- phyf::bird_beak_codes |> pf_as_phylo()

# For each internal node, compute mean curvature of all descendant edges
n_tips <- length(phylo_tree$tip.label)
n_nodes <- phylo_tree$Nnode

# For each internal node, find descendant tips via ape::extract.clade
# and compute mean curvature of those tip edges
get_desc_tips <- function(node_num, tree) {
  subtree <- ape::extract.clade(tree, node_num)
  subtree$tip.label
}

clade_scores <- tibble(
  node_num = (n_tips + 1):(n_tips + n_nodes),
  node_label = phylo_tree$node.label
) |>
  mutate(
    desc_tips = map(node_num,
                    ~ get_desc_tips(.x, phylo_tree)),
    n_desc_tips = map_int(desc_tips, length)
  ) |>
  # Skip root and very large clades
  filter(n_desc_tips >= 2, n_desc_tips <= 500) |>
  mutate(
    mean_curvature = map_dbl(desc_tips, function(dt) {
      desc_scores <- edge_scores |>
        filter(edge %in% dt)
      if (nrow(desc_scores) == 0) return(NA_real_)
      mean(desc_scores$curvature_ratio, na.rm = TRUE)
    }),
    max_curvature = map_dbl(desc_tips, function(dt) {
      desc_scores <- edge_scores |>
        filter(edge %in% dt)
      if (nrow(desc_scores) == 0) return(NA_real_)
      max(desc_scores$curvature_ratio, na.rm = TRUE)
    })
  ) |>
  arrange(desc(mean_curvature))

message("\n=== TOP 15 CLADES BY MEAN TIP CURVATURE ===")
print(
  clade_scores |>
    select(node_label, n_desc_tips, mean_curvature, max_curvature) |>
    head(15),
  n = 15
)

write_csv(
  clade_scores |> select(-desc_tips),
  "output/edge_scan/clade_curvature_scores.csv"
)

message("\nDone! Results saved to output/edge_scan/")
