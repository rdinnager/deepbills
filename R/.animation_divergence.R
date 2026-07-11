## Phase 4+5: Divergence Animation & Manifold vs Linear Comparison
## GPU NEEDED: Yes (SDF model + VAE decoder)
##
## Two species beaks side-by-side, animated from tips toward their
## common ancestor (beaks converge), then reversed (beaks diverge).
## Also generates manifold vs straightened comparison for selected edges.
##
## Run AFTER .animation_edge_scan.R for species selection.

options(rgl.useNULL = TRUE)
options(torch.serialization_version = 2)

library(tidyverse)
library(torch)
library(fibre)
library(phyf)
library(ape)
library(rgl)
library(magick)

source("R/.animation_utils.R")

# ---- Load data ----

message("Loading data...")
z_mani <- read_rds("data/evo_predictions_all_edges_v2.rds")
scaling <- read_rds("data/code_dat_means_sds.rds")
active_dims <- read_rds("data/active_dims_16dim.rds")
niche_pal <- load_niche_pal()
trophic_dat <- read_rds("data/trophic_dat.rds")
bill_edge_trajs <- read_rds("data/bill_edge_trajs_16dim.rds")

# Build phylo tree
phylo_tree <- phyf::bird_beak_codes |> pf_as_phylo()

# Load models
message("Loading SDF model...")
beak_mod <- load_bird_beak_model()
beak_mod <- beak_mod$cuda()

# ---- Path finding utilities ----

#' Map node numbers to edge labels in z_mani
#' @param node_num integer node number from phylo_tree
#' @param tree phylo object
node_to_label <- function(node_num, tree) {
  n_tips <- length(tree$tip.label)
  if (node_num <= n_tips) {
    tree$tip.label[node_num]
  } else {
    tree$node.label[node_num - n_tips]
  }
}

#' Get the sequence of edge labels from a tip to the MRCA
#' @param tip_name species name
#' @param mrca_node node number of MRCA
#' @param tree phylo object
#' @return character vector of edge labels (tip first, MRCA last)
get_path_edges <- function(tip_name, mrca_node, tree) {
  tip_idx <- which(tree$tip.label == tip_name)
  path_nodes <- ape::nodepath(tree, tip_idx, mrca_node)
  # Convert node numbers to labels
  map_chr(path_nodes, ~ node_to_label(.x, tree))
}

#' Extract concatenated trajectory along a path (tip → MRCA)
#' Each edge's z_seqs/prcodes/prtroph are reversed since they
#' are stored parent→child but we traverse child→parent.
#' @param edge_labels ordered vector of edge labels (tip first)
#' @param z_tree_df the predictions tibble
#' @param subsample_every take every Nth point (default 1 = all)
get_path_trajectory <- function(edge_labels, z_tree_df,
                                 subsample_every = 1) {
  # Collect data for each edge in path order
  path_data <- map(edge_labels, function(lab) {
    row <- z_tree_df |> filter(edge == lab)
    if (nrow(row) == 0) return(NULL)
    list(
      z_seqs = row$z_seqs[[1]],
      prcodes = row$prcodes[[1]],
      prtroph = row$prtroph[[1]]
    )
  }) |>
    compact()

  # Reverse each edge (go from child/tip toward parent/root)
  # and concatenate
  all_prcodes <- map(path_data, ~ {
    .x$prcodes[nrow(.x$prcodes):1, , drop = FALSE]
  }) |> do.call(rbind, args = _)

  all_prtroph <- map(path_data, ~ rev(.x$prtroph)) |> unlist()

  all_z_seqs <- map(path_data, ~ {
    .x$z_seqs[nrow(.x$z_seqs):1, ]
  }) |> bind_rows()

  # Subsample
  if (subsample_every > 1) {
    idx <- seq(1, nrow(all_prcodes), by = subsample_every)
    all_prcodes <- all_prcodes[idx, , drop = FALSE]
    all_prtroph <- all_prtroph[idx]
    all_z_seqs <- all_z_seqs[idx, ]
  }

  list(
    prcodes = all_prcodes,
    prtroph = all_prtroph,
    z_seqs = all_z_seqs,
    n_frames = nrow(all_prcodes)
  )
}

# ---- Species selection ----

# Set species pair here (or choose from edge_scores)
# Examples - uncomment or set your own:
# species1 <- "Archilochus_colubris"   # Ruby-throated Hummingbird
# species2 <- "Pelecanus_occidentalis" # Brown Pelican

# Default: pick top 2 tip edges by curvature ratio
edge_scores <- read_csv(
  "output/edge_scan/edge_curvature_scores.csv",
  show_col_types = FALSE
)

top_tips <- edge_scores |>
  filter(is_tip, curvature_ratio > 1) |>
  slice_max(curvature_ratio, n = 10)

message("\nTop 10 tip edges by curvature ratio:")
print(top_tips |> select(edge, curvature_ratio, transition),
      n = 10)

# Pick first two by default
if (!exists("species1")) {
  species1 <- top_tips$edge[1]
}
if (!exists("species2")) {
  species2 <- top_tips$edge[2]
}

message("\nSelected species pair:")
message("  Species 1: ", species1)
message("  Species 2: ", species2)

# ---- Find MRCA and paths ----

mrca_node <- ape::getMRCA(
  phylo_tree,
  c(which(phylo_tree$tip.label == species1),
    which(phylo_tree$tip.label == species2))
)
mrca_label <- node_to_label(mrca_node, phylo_tree)

message("MRCA: ", mrca_label, " (node ", mrca_node, ")")

path1_edges <- get_path_edges(species1, mrca_node, phylo_tree)
path2_edges <- get_path_edges(species2, mrca_node, phylo_tree)

message("Path 1 (", species1, " -> MRCA): ",
        length(path1_edges), " nodes")
message("Path 2 (", species2, " -> MRCA): ",
        length(path2_edges), " nodes")

# ---- Extract trajectories ----

# Subsample to keep animation manageable
# Target ~100-150 frames per path
path1_total <- length(path1_edges) * 50
path2_total <- length(path2_edges) * 50
subsample1 <- max(1, path1_total %/% 120)
subsample2 <- max(1, path2_total %/% 120)

traj1 <- get_path_trajectory(path1_edges, z_mani, subsample1)
traj2 <- get_path_trajectory(path2_edges, z_mani, subsample2)

message("Path 1 frames: ", traj1$n_frames,
        " (subsampled every ", subsample1, ")")
message("Path 2 frames: ", traj2$n_frames,
        " (subsampled every ", subsample2, ")")

# ---- Equalize frame counts ----

# Both paths should have the same number of frames
# for synchronized side-by-side animation
n_sync <- max(traj1$n_frames, traj2$n_frames)

resample_trajectory <- function(traj, n_target) {
  if (traj$n_frames == n_target) return(traj)
  idx <- round(seq(1, traj$n_frames, length.out = n_target))
  list(
    prcodes = traj$prcodes[idx, , drop = FALSE],
    prtroph = traj$prtroph[idx],
    z_seqs = traj$z_seqs[idx, ],
    n_frames = n_target
  )
}

traj1 <- resample_trajectory(traj1, n_sync)
traj2 <- resample_trajectory(traj2, n_sync)

message("Synchronized to ", n_sync, " frames each")

# ---- Generate meshes for both paths ----

message("\nGenerating meshes for path 1 (", species1, ")...")
meshes1 <- map(seq_len(n_sync), function(i) {
  code <- matrix(traj1$prcodes[i, ], nrow = 1)
  make_mesh(code, scaling, beak_mod, resolution = 100)
}, .progress = TRUE)

message("Generating meshes for path 2 (", species2, ")...")
meshes2 <- map(seq_len(n_sync), function(i) {
  code <- matrix(traj2$prcodes[i, ], nrow = 1)
  make_mesh(code, scaling, beak_mod, resolution = 100)
}, .progress = TRUE)

# ---- Compute colors ----

colors1 <- compute_smooth_colors(traj1$prtroph, niche_pal, 5)
colors2 <- compute_smooth_colors(traj2$prtroph, niche_pal, 5)

# ---- Render divergence animation ----

outdir <- "output/animations"
frame_files <- character(n_sync)

message("\nRendering ", n_sync, " side-by-side frames...")
for (i in seq_len(n_sync)) {
  # Render left beak
  left_file <- tempfile(fileext = ".png")
  render_mesh_frame(meshes1[[i]], col = colors1[i],
                     filename = left_file,
                     width = 400, height = 400)

  # Render right beak
  right_file <- tempfile(fileext = ".png")
  render_mesh_frame(meshes2[[i]], col = colors2[i],
                     filename = right_file,
                     width = 400, height = 400)

  # Combine side by side
  left_img <- image_read(left_file)
  right_img <- image_read(right_file)
  combined <- image_append(c(left_img, right_img))

  # Add species labels
  sp1_display <- gsub("_", " ", species1)
  sp2_display <- gsub("_", " ", species2)

  combined <- image_annotate(
    combined, sp1_display,
    size = 20, color = "black",
    gravity = "northwest", font = "Helvetica",
    location = "+10+10"
  )
  combined <- image_annotate(
    combined, sp2_display,
    size = 20, color = "black",
    gravity = "northeast", font = "Helvetica",
    location = "+10+10"
  )

  # Trophic niche labels
  combined <- image_annotate(
    combined, traj1$prtroph[i],
    size = 18, color = colors1[i],
    gravity = "southwest", font = "Helvetica-Bold",
    location = "+10+10",
    strokecolor = "white", strokewidth = 1
  )
  combined <- image_annotate(
    combined, traj2$prtroph[i],
    size = 18, color = colors2[i],
    gravity = "southeast", font = "Helvetica-Bold",
    location = "+10+10",
    strokecolor = "white", strokewidth = 1
  )

  # Progress indicator
  progress_pct <- round(100 * i / n_sync)
  combined <- image_annotate(
    combined, paste0("-> MRCA (", progress_pct, "%)"),
    size = 14, color = "grey40",
    gravity = "south", font = "Helvetica",
    location = "+0+35"
  )

  fname <- file.path(outdir, sprintf("diverge_%03d.png", i))
  image_write(combined, fname)
  frame_files[i] <- fname

  # Cleanup temp files
  unlink(c(left_file, right_file))

  if (i %% 20 == 0) message("  Frame ", i, "/", n_sync)
}

# ---- Assemble GIF: convergence + hold + divergence ----

hold_frames <- rep(frame_files[n_sync], 20)  # hold at MRCA
convergence <- frame_files
divergence <- rev(frame_files[-c(1, n_sync)])

all_frames <- c(convergence, hold_frames, divergence)

gif_name <- sprintf("divergence_%s_vs_%s.gif",
                     substr(species1, 1, 20),
                     substr(species2, 1, 20))
gif_file <- file.path(outdir, gif_name)

assemble_gif(all_frames, gif_file, fps = 15,
             width = 800, height = 400)

message("\nDivergence animation complete!")
message("GIF: ", gif_file)

# ============================================================
# Phase 5: Manifold vs Straightened Comparison
# ============================================================

message("\n=== PHASE 5: Manifold vs Straightened Comparison ===")

# Pick the highest-curvature edge from the pair's paths
path_edges_all <- unique(c(path1_edges, path2_edges))
comparison_edge <- edge_scores |>
  filter(edge %in% path_edges_all) |>
  slice_max(curvature_ratio, n = 1) |>
  pull(edge)

message("Comparison edge: ", comparison_edge,
        " (curvature ratio: ",
        round(edge_scores$curvature_ratio[
          edge_scores$edge == comparison_edge
        ], 2), ")")

# ---- Get manifold trajectory for this edge ----

edge_row <- z_mani |> filter(edge == comparison_edge)
manifold_prcodes <- edge_row$prcodes[[1]]
manifold_prtroph <- edge_row$prtroph[[1]]

# ---- Generate straightened trajectory ----
# Linear interpolation between start and end z_seqs,
# then decode through VAE

message("Loading VAE for straightened path decoding...")
vae <- torch_load("data/bill_vae_w_trophic_v1.to")
vae <- vae$cuda()

z_start <- as.matrix(edge_row$z_seqs[[1]])[1, ]
z_end <- as.matrix(edge_row$z_seqs[[1]])[50, ]

# Linearly interpolate 50 points
t_vals <- seq(0, 1, length.out = 50)
straight_z_seqs <- t(sapply(t_vals, function(t) {
  (1 - t) * z_start + t * z_end
}))

# Decode through VAE
trophic_levs <- levels(trophic_dat$trophic_niche)

decode_zseqs <- function(zseq_mat, active_dims, trophic_levs) {
  z_mat <- matrix(0, nrow = nrow(zseq_mat), ncol = 64)
  z_mat[, active_dims] <- zseq_mat
  with_no_grad({
    z_tensor <- torch_tensor(z_mat, device = "cuda")
    x <- vae$decoder(z_tensor)
    pred_troph <- torch_argmax(
      nnf_softmax(x$out_trophic, 2L), dim = 2L
    )
    pred_codes <- x$out_codes
  })
  list(
    troph = trophic_levs[as.matrix(pred_troph$cpu())],
    codes = as.matrix(pred_codes$cpu())
  )
}

message("Decoding straightened path...")
straight_decoded <- decode_zseqs(
  straight_z_seqs, active_dims, trophic_levs
)
straight_prcodes <- straight_decoded$codes
straight_prtroph <- straight_decoded$troph

# ---- Generate meshes for both versions ----

frame_idx <- seq(1, 50, by = 2)  # 25 frames each
n_comp <- length(frame_idx)

message("Generating manifold meshes...")
meshes_mani <- map(frame_idx, function(i) {
  code <- matrix(manifold_prcodes[i, ], nrow = 1)
  make_mesh(code, scaling, beak_mod, resolution = 100)
}, .progress = TRUE)

message("Generating straightened meshes...")
meshes_straight <- map(frame_idx, function(i) {
  code <- matrix(straight_prcodes[i, ], nrow = 1)
  make_mesh(code, scaling, beak_mod, resolution = 100)
}, .progress = TRUE)

# Colors
colors_mani <- compute_smooth_colors(
  manifold_prtroph[frame_idx], niche_pal, 3
)
colors_straight <- compute_smooth_colors(
  straight_prtroph[frame_idx], niche_pal, 3
)

# ---- Render comparison frames ----

comp_files <- character(n_comp)

message("Rendering ", n_comp, " comparison frames...")
for (i in seq_len(n_comp)) {
  # Top: manifold
  top_file <- tempfile(fileext = ".png")
  render_mesh_frame(meshes_mani[[i]], col = colors_mani[i],
                     filename = top_file, width = 600, height = 400)

  # Bottom: straightened
  bot_file <- tempfile(fileext = ".png")
  render_mesh_frame(meshes_straight[[i]], col = colors_straight[i],
                     filename = bot_file, width = 600, height = 400)

  top_img <- image_read(top_file)
  bot_img <- image_read(bot_file)

  # Labels
  top_img <- image_annotate(
    top_img, "Manifold (curved path)",
    size = 20, color = "black",
    gravity = "northwest", font = "Helvetica",
    location = "+10+10"
  )
  top_img <- image_annotate(
    top_img, manifold_prtroph[frame_idx[i]],
    size = 16, color = colors_mani[i],
    gravity = "southwest", font = "Helvetica-Bold",
    location = "+10+10",
    strokecolor = "white", strokewidth = 1
  )

  bot_img <- image_annotate(
    bot_img, "Straightened (linear interpolation)",
    size = 20, color = "black",
    gravity = "northwest", font = "Helvetica",
    location = "+10+10"
  )
  bot_img <- image_annotate(
    bot_img, straight_prtroph[frame_idx[i]],
    size = 16, color = colors_straight[i],
    gravity = "southwest", font = "Helvetica-Bold",
    location = "+10+10",
    strokecolor = "white", strokewidth = 1
  )

  combined <- image_append(c(top_img, bot_img), stack = TRUE)

  fname <- file.path(outdir,
                      sprintf("compare_%s_%03d.png",
                              comparison_edge, i))
  image_write(combined, fname)
  comp_files[i] <- fname

  unlink(c(top_file, bot_file))
}

# Boomerang GIF
comp_all <- c(comp_files, rev(comp_files[-c(1, n_comp)]))
comp_gif <- file.path(outdir,
                       sprintf("comparison_%s.gif", comparison_edge))
assemble_gif(comp_all, comp_gif, fps = 12,
             width = 600, height = 800)

message("\nComparison animation complete!")
message("GIF: ", comp_gif)
message("Edge: ", comparison_edge)

message("\n=== ALL DONE ===")
message("Divergence GIF: ", gif_file)
message("Comparison GIF: ", comp_gif)
