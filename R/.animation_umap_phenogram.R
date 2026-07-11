## Phase 3: UMAP Phenogram — 2D latent space tree with beak images
## GPU NEEDED: Yes (for rendering beak thumbnails at selected tips)
##
## Projects all evolutionary trajectories into 2D UMAP space,
## draws the phylogenetic tree as colored paths, and overlays
## small rendered beak images at selected tip positions.
## Uses true observed latent codes for tip beaks (not model estimates).

options(rgl.useNULL = TRUE)
options(torch.serialization_version = 2)

library(tidyverse)
library(uwot)
library(torch)
library(fibre)
library(rgl)
library(ape)
library(phyf)
library(targets)

source("R/.animation_utils.R")

# ---- Load data ----

message("Loading data...")
z_mani <- read_rds("data/evo_predictions_all_edges_v2.rds")
niche_pal <- load_niche_pal()
trophic_dat <- read_rds("data/trophic_dat.rds")
avonet <- tar_read(bird_beak_avonet)
scaling <- read_rds("data/code_dat_means_sds.rds")

# ---- 3A. UMAP embedding ----

umap_file <- "output/edge_scan/umap_2d_embedding.rds"

if (file.exists(umap_file)) {
  message("Loading existing UMAP embedding...")
  umap_result <- read_rds(umap_file)
  embedding_2d <- umap_result$embedding
} else {
  message("Stacking z_seqs for UMAP (", nrow(z_mani),
          " edges x 50 points)...")
  big_z_mat <- do.call(rbind, lapply(z_mani$z_seqs, as.matrix))
  message("Matrix size: ", nrow(big_z_mat), " x ", ncol(big_z_mat))

  message("Running UMAP (n_components=2)...")
  embedding_2d <- uwot::umap(big_z_mat, n_components = 2,
                              n_neighbors = 30, min_dist = 0.3)
  umap_result <- list(
    embedding = embedding_2d,
    edges = z_mani$edge,
    n_points_per_edge = 50
  )
  write_rds(umap_result, umap_file)
  message("UMAP embedding saved.")
}

# ---- 3B. Build per-edge trajectory dataframe ----

umap_df <- tibble(
  edge = rep(z_mani$edge, each = 50),
  is_tip = rep(z_mani$is_tip, each = 50),
  prtroph = unlist(z_mani$prtroph),
  U1 = embedding_2d[, 1],
  U2 = embedding_2d[, 2],
  time_point = rep(1:50, nrow(z_mani))
)

tip_endpoints <- umap_df |>
  filter(is_tip, time_point == 50) |>
  distinct(edge, .keep_all = TRUE) |>
  left_join(
    trophic_dat |> select(edge = label, true_trophic = trophic_niche),
    by = "edge"
  ) |>
  left_join(
    avonet |> select(edge = label, Order, BLFamilyEnglish),
    by = "edge"
  )

# ---- Order common names ----

order_common <- tribble(
  ~Order, ~clade_label,
  "ACCIPITRIFORMES", "Birds of Prey",
  "ANSERIFORMES", "Ducks &\nGeese",
  "APODIFORMES", "Hummingbirds\n& Swifts",
  "CHARADRIIFORMES", "Shorebirds",
  "COLUMBIFORMES", "Pigeons &\nDoves",
  "CORACIIFORMES", "Kingfishers",
  "GALLIFORMES", "Gamebirds",
  "GRUIFORMES", "Rails &\nCranes",
  "PELECANIFORMES", "Pelicans &\nHerons",
  "PICIFORMES", "Woodpeckers\n& Toucans",
  "PSITTACIFORMES", "Parrots",
  "STRIGIFORMES", "Owls",
  "PROCELLARIIFORMES", "Albatrosses\n& Petrels",
  "SPHENISCIFORMES", "Penguins",
  "FALCONIFORMES", "Falcons",
  "BUCEROTIFORMES", "Hornbills",
  "CUCULIFORMES", "Cuckoos"
)

# ---- 3C. Select tips: trophic-aware + spatial coverage ----

message("Selecting tips for annotation...")

# Convex hull tips for edge coverage
hull_idx <- chull(tip_endpoints$U1, tip_endpoints$U2)
hull_tips <- tip_endpoints$edge[hull_idx]

# Per-trophic-niche representatives in interior
# Pick well-spaced tips within each niche cluster
set.seed(42)
interior <- tip_endpoints |> filter(!edge %in% hull_tips)

niche_reps <- interior |>
  group_by(true_trophic) |>
  group_modify(function(df, key) {
    n_pick <- min(7, nrow(df))
    if (n_pick <= 1) return(df)
    km <- kmeans(cbind(df$U1, df$U2),
                 centers = n_pick, nstart = 5)
    df |>
      mutate(cluster = km$cluster) |>
      group_by(cluster) |>
      slice_sample(n = 1) |>
      ungroup()
  }) |>
  ungroup()

# Extra picks from specific clades of interest
extra_clades <- c("APODIFORMES", "PICIFORMES")
extra_picks <- tip_endpoints |>
  filter(Order %in% extra_clades,
         !edge %in% c(hull_tips, niche_reps$edge)) |>
  group_by(Order) |>
  group_modify(function(df, key) {
    n_pick <- min(6, nrow(df))
    if (n_pick <= 1) return(df)
    km <- kmeans(cbind(df$U1, df$U2),
                 centers = n_pick, nstart = 5)
    df |>
      mutate(cluster = km$cluster) |>
      group_by(cluster) |>
      slice_sample(n = 1) |>
      ungroup()
  }) |>
  ungroup()

# Additional interior fill: k-means on ALL remaining tips
already_picked <- c(hull_tips, niche_reps$edge, extra_picks$edge)
remaining <- tip_endpoints |>
  filter(!edge %in% already_picked)
set.seed(123)
fill_km <- kmeans(cbind(remaining$U1, remaining$U2),
                  centers = 20, nstart = 5)
fill_picks <- remaining |>
  mutate(cluster = fill_km$cluster) |>
  group_by(cluster) |>
  slice_sample(n = 1) |>
  ungroup()

annotated_tips <- unique(c(
  hull_tips, niche_reps$edge,
  extra_picks$edge, fill_picks$edge
))
message("  Hull: ", length(hull_tips),
        ", Niche reps: ", nrow(niche_reps),
        ", Total: ", length(annotated_tips))

write_rds(annotated_tips, "output/edge_scan/umap_annotated_tips.rds")

# ---- 3D. Render thumbnails with TRUE latent codes ----

message("Loading SDF model...")
beak_mod <- load_bird_beak_model()
beak_mod <- beak_mod$cuda()

thumb_dir <- "output/animations/tip_thumbs"

# Get true 64-dim latent codes from AVONET (not model estimates)
true_codes <- avonet |>
  select(edge = label, starts_with("latent_code_")) |>
  filter(edge %in% annotated_tips)

message("Rendering ", nrow(true_codes),
        " thumbnails with true latent codes...")
for (i in seq_len(nrow(true_codes))) {
  tip_edge <- true_codes$edge[i]
  code_vec <- true_codes |>
    filter(edge == tip_edge) |>
    select(starts_with("latent_code_")) |>
    as.matrix()

  troph <- tip_endpoints$true_trophic[
    tip_endpoints$edge == tip_edge
  ]
  if (length(troph) == 0 || is.na(troph)) {
    troph <- "Invertivore"
  }

  fname <- file.path(thumb_dir, paste0(tip_edge, ".png"))

  mesh <- tryCatch(
    beak_mod$get_mesh(
      torch_tensor(code_vec, device = "cuda"),
      resolution = 80, smooth = FALSE
    ),
    error = function(e) NULL
  )

  if (!is.null(mesh)) {
    render_mesh_frame(
      mesh, col = niche_pal[troph],
      filename = fname,
      width = 200, height = 200,
      view = "3/4",
      transparent_bg = TRUE
    )
  }

  if (i %% 20 == 0) {
    message("  ", i, "/", nrow(true_codes))
  }
}
message("Thumbnails done!")

# ---- 3E. Clade label positions ----

# Compute centroid of each order's tips in UMAP space
clade_labels <- tip_endpoints |>
  inner_join(order_common, by = "Order") |>
  group_by(clade_label) |>
  summarise(
    U1 = median(U1),
    U2 = median(U2),
    n = n(),
    .groups = "drop"
  ) |>
  filter(n >= 5)

# Use ggrepel to compute non-overlapping label positions.
# Beak thumbnail positions are included as empty-label obstacles
# so that real clade labels are repelled away from them.
tip_pos <- tip_endpoints |>
  filter(edge %in% annotated_tips) |>
  select(U1, U2)

# Find non-overlapping label positions by searching for the nearest
# clear spot around each clade centroid, avoiding beak bounding boxes
# and previously placed labels.

tip_pos <- tip_endpoints |>
  filter(edge %in% annotated_tips) |>
  select(U1, U2)

# Beak bounding boxes in data coordinates
u1_range <- diff(range(umap_df$U1, na.rm = TRUE))
u2_range <- diff(range(umap_df$U2, na.rm = TRUE))
offset <- min(u1_range, u2_range) * 0.028
beak_half <- offset

# Check if a candidate position overlaps any beak or placed label
overlaps_any <- function(cx, cy, half_w, half_h,
                          boxes) {
  if (nrow(boxes) == 0) return(FALSE)
  any(
    cx - half_w < boxes$xmax &
    cx + half_w > boxes$xmin &
    cy - half_h < boxes$ymax &
    cy + half_h > boxes$ymin
  )
}

# Build obstacle list: beak thumbnail bounding boxes
obstacles <- tibble(
  xmin = tip_pos$U1 - beak_half,
  xmax = tip_pos$U1 + beak_half,
  ymin = tip_pos$U2 - beak_half,
  ymax = tip_pos$U2 + beak_half
)

# Approximate label size in data units
label_half_w <- 1.2
label_half_h <- 0.6

# For each clade label, spiral outward from centroid until clear
clade_labels$label_U1 <- NA_real_
clade_labels$label_U2 <- NA_real_

for (i in seq_len(nrow(clade_labels))) {
  cx <- clade_labels$U1[i]
  cy <- clade_labels$U2[i]
  placed <- FALSE

  # Try positions in expanding circles
  for (radius in seq(0, 6, by = 0.3)) {
    if (radius == 0) {
      angles <- 0
    } else {
      angles <- seq(0, 2 * pi, length.out = max(8, round(radius * 4)))
    }
    for (a in angles) {
      test_x <- cx + radius * cos(a)
      test_y <- cy + radius * sin(a)

      if (!overlaps_any(test_x, test_y,
                        label_half_w, label_half_h,
                        obstacles)) {
        clade_labels$label_U1[i] <- test_x
        clade_labels$label_U2[i] <- test_y

        # Add this label as an obstacle for subsequent labels
        obstacles <- bind_rows(obstacles, tibble(
          xmin = test_x - label_half_w,
          xmax = test_x + label_half_w,
          ymin = test_y - label_half_h,
          ymax = test_y + label_half_h
        ))
        placed <- TRUE
        break
      }
    }
    if (placed) break
  }

  # Fallback: use centroid
  if (!placed) {
    clade_labels$label_U1[i] <- cx
    clade_labels$label_U2[i] <- cy
  }
}

message("Clade labels: ", nrow(clade_labels),
        " (positions computed via ggrepel)")

# ---- 3F. Build phenogram plot ----

message("Building plot...")

tip_img_df <- tip_endpoints |>
  filter(edge %in% annotated_tips) |>
  mutate(image = file.path(thumb_dir, paste0(edge, ".png"))) |>
  filter(file.exists(image))

message("  ", nrow(tip_img_df), " thumbnails available")

p <- ggplot() +
  # Edge trajectories
  geom_path(
    data = umap_df,
    aes(x = U1, y = U2, group = edge, color = prtroph),
    alpha = 0.25, linewidth = 0.3
  ) +
  # Tip points
  geom_point(
    data = tip_endpoints,
    aes(x = U1, y = U2, color = true_trophic),
    size = 0.5, alpha = 1
  ) +
  scale_color_manual(values = niche_pal, name = "Trophic Niche") +
  # Clade labels at ggrepel-computed positions
  # Leader lines from centroid to label
  geom_segment(
    data = clade_labels,
    aes(x = U1, y = U2, xend = label_U1, yend = label_U2),
    color = "grey50", linewidth = 0.3, alpha = 0.6
  ) +
  geom_label(
    data = clade_labels,
    aes(x = label_U1, y = label_U2, label = clade_label),
    size = 3, alpha = 0.8, fill = "white",
    label.padding = unit(0.15, "lines"),
    linewidth = 0.2,
    fontface = "bold"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14, face = "bold"),
    legend.key.size = unit(1.2, "lines"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 16)
  ) +
  guides(
    color = guide_legend(
      override.aes = list(size = 4, alpha = 1)
    )
  ) +
  labs(
    x = "UMAP 1", y = "UMAP 2",
    title = "Beak Shape Evolution on the Bird Phylogeny",
    subtitle = "Edges colored by predicted trophic niche"
  )

# Add beak images
for (i in seq_len(nrow(tip_img_df))) {
  img <- png::readPNG(tip_img_df$image[i])
  p <- p +
    annotation_custom(
      grid::rasterGrob(img, interpolate = TRUE),
      xmin = tip_img_df$U1[i] - offset,
      xmax = tip_img_df$U1[i] + offset,
      ymin = tip_img_df$U2[i] - offset,
      ymax = tip_img_df$U2[i] + offset
    )
}

# Save base (no beaks) and full versions
ggsave("output/animations/umap_phenogram_base.pdf",
       p %+% list(), width = 16, height = 13)

ggsave("output/animations/umap_phenogram_with_beaks.pdf",
       p, width = 16, height = 13)
ggsave("output/animations/umap_phenogram_with_beaks.png",
       p, width = 16, height = 13, dpi = 200)

message("\nPhase 3 complete!")
message("Outputs in output/animations/umap_phenogram_*")
