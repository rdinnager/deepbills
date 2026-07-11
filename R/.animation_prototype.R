## Phase 2: Prototype — Single Edge Beak Morph GIF
## GPU NEEDED: Yes (SDF model + VAE for linear paths)
##
## Generates animated GIFs for selected high-curvature edges,
## showing the beak morphing along each branch, colored by trophic niche.
## Produces both manifold (curved) and linear (straightened) animations.
##
## Run AFTER .animation_edge_scan.R to have edge scores available.

options(rgl.useNULL = TRUE)
options(torch.serialization_version = 2)

library(tidyverse)
library(torch)
library(fibre)
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
trophic_levs <- levels(trophic_dat$trophic_niche)

edge_scores <- read_csv(
  "output/edge_scan/edge_curvature_scores.csv",
  show_col_types = FALSE
)

# ---- Select edges ----

test_edges <- c(
  "Phainoptila_melanoxantha",   # Invertivore -> Frugivore
  "Myadestes_occidentalis",     # Invertivore -> Frugivore
  "Rostratula_benghalensis",    # Omnivore -> Aquatic predator
  "Leptopoecile_sophiae"        # No transition, highest excess
)

message("Test edges:")
edge_scores |>
  filter(edge %in% test_edges) |>
  select(edge, curvature_ratio, excess_dist,
         has_transition, transition) |>
  print(n = 4, width = 100)

# ---- Load models ----

message("\nLoading SDF model...")
beak_mod <- load_bird_beak_model()
beak_mod <- beak_mod$cuda()

message("Loading VAE model...")
vae <- torch_load("data/bill_vae_w_trophic_v1.to")
vae <- vae$cuda()

# ---- VAE decoder ----

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

# ---- Core animation function ----

#' Generate a prototype animation GIF for one edge
#' @param mode "manifold" uses the curved path from the model,
#'             "linear" linearly interpolates the latent endpoints
make_prototype <- function(chosen_edge, z_mani, scaling,
                            beak_mod, niche_pal,
                            edge_scores,
                            mode = "manifold",
                            outdir = "output/animations",
                            resolution = 100,
                            fps = 15) {
  label <- ifelse(mode == "linear", "LINEAR", "MANIFOLD")
  message("\n========================================")
  message(label, ": ", chosen_edge)

  score_row <- edge_scores |>
    filter(edge == chosen_edge)
  message("  Curvature ratio: ",
          round(score_row$curvature_ratio, 3))
  message("  Excess distance: ",
          round(score_row$excess_dist, 3))
  if (!is.na(score_row$transition)) {
    message("  Transition: ", score_row$transition)
  }

  # Extract edge data
  edge_row <- z_mani |> filter(edge == chosen_edge)

  if (mode == "manifold") {
    # Use pre-computed curved path
    prcodes <- edge_row$prcodes[[1]]
    prtroph <- edge_row$prtroph[[1]]
  } else {
    # Linear interpolation in 16-dim latent space
    z_seqs <- as.matrix(edge_row$z_seqs[[1]])
    z_start <- z_seqs[1, ]
    z_end <- z_seqs[nrow(z_seqs), ]
    t_vals <- seq(0, 1, length.out = 50)
    z_linear <- t(sapply(t_vals, function(t) {
      (1 - t) * z_start + t * z_end
    }))

    message("  Decoding linearized path through VAE...")
    decoded <- decode_zseqs(z_linear, active_dims, trophic_levs)
    prcodes <- decoded$codes
    prtroph <- decoded$troph
  }

  message("  Trophic: ",
          paste(unique(prtroph), collapse = " -> "))

  # Arc-length resample for smooth animation
  resampled <- resample_by_arclength(
    prcodes, n_out = nrow(prcodes), troph = prtroph
  )
  prcodes <- resampled$codes
  prtroph <- resampled$troph

  n_frames <- nrow(prcodes)

  message("  Generating ", n_frames,
          " meshes at resolution ", resolution, "...")
  meshes <- map(seq_len(n_frames), function(i) {
    code <- matrix(prcodes[i, ], nrow = 1)
    make_mesh(code, scaling, beak_mod,
              resolution = resolution)
  }, .progress = TRUE)

  # Smoothed trophic colors
  frame_colors <- compute_smooth_colors(
    prtroph, niche_pal, blend_window = 7
  )

  # Render frames
  prefix <- ifelse(mode == "linear", "proto_lin", "proto")
  frame_files <- character(n_frames)
  message("  Rendering ", n_frames, " frames...")
  for (i in seq_len(n_frames)) {
    fname <- file.path(
      outdir,
      sprintf("%s_%s_%03d.png", prefix, chosen_edge, i)
    )
    render_mesh_frame(
      meshes[[i]], col = frame_colors[i],
      filename = fname,
      width = 800, height = 600,
      view = "3/4"
    )
    frame_files[i] <- fname
  }

  # Labels + timeline
  mode_label <- ifelse(mode == "linear", " (linear)", "")
  species_display <- paste0(
    gsub("_", " ", chosen_edge), mode_label
  )
  message("  Adding labels + timeline...")
  for (i in seq_len(n_frames)) {
    img <- image_read(frame_files[i])

    img <- image_annotate(
      img, species_display,
      size = 24, color = "black",
      gravity = "north", font = "Helvetica",
      location = "+0+10"
    )

    img <- image_annotate(
      img, prtroph[i],
      size = 22, color = frame_colors[i],
      gravity = "south", font = "Helvetica-Bold",
      location = "+0+80",
      strokecolor = "white", strokewidth = 1
    )

    progress <- (i - 1) / (n_frames - 1)
    img <- add_timeline(
      img, progress, width = 800,
      y_pos = 40, bar_color = "grey40",
      marker_color = frame_colors[i]
    )

    image_write(img, frame_files[i])
  }

  # GIF with boomerang
  suffix <- ifelse(mode == "linear", "_linear", "")
  gif_file <- file.path(
    outdir,
    sprintf("prototype%s_%s.gif", suffix, chosen_edge)
  )
  boomerang <- c(
    frame_files,
    rev(frame_files[-c(1, n_frames)])
  )
  assemble_gif(boomerang, gif_file,
               fps = fps, width = 800, height = 600)

  # Clean up frame PNGs
  unlink(frame_files)
  message("  Cleaned up ", length(frame_files), " frame PNGs")

  list(
    edge = chosen_edge,
    mode = mode,
    gif = gif_file,
    n_frames = n_frames,
    transition = paste(unique(prtroph), collapse = " -> ")
  )
}

# ---- Run both manifold and linear for all test edges ----

results <- list()
for (e in test_edges) {
  for (m in c("manifold", "linear")) {
    r <- tryCatch(
      make_prototype(
        e, z_mani, scaling, beak_mod,
        niche_pal, edge_scores, mode = m
      ),
      error = function(err) {
        message("  ERROR: ", err$message)
        list(edge = e, mode = m, gif = NA,
             error = err$message)
      }
    )
    results <- c(results, list(r))
  }
}

# ---- Summary ----

message("\n========================================")
message("=== PROTOTYPE RESULTS ===")
message("========================================")
for (r in results) {
  if (!is.null(r$error)) {
    message("FAILED: ", r$edge, " (", r$mode, ") - ",
            r$error)
  } else {
    message(
      r$edge, " [", r$mode, "]",
      "  trophic=", r$transition
    )
    message("  -> ", r$gif)
  }
}

message("\nDone!")
