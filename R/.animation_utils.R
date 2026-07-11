## Shared utilities for beak evolution animations
## Used by: .animation_edge_scan.R, .animation_prototype.R,
##          .animation_umap_phenogram.R, .animation_divergence.R

# IMPORTANT: rgl.useNULL must be set BEFORE library(rgl)
# to avoid segfaults on headless Linux
options(rgl.useNULL = TRUE)

# Chrome is needed for rgl snapshot rendering via webshot2
# On HiPerGator: module load chrome -> /apps/ubuntu/24.04/bin/chrome
if (Sys.getenv("CHROMOTE_CHROME") == "" &&
    file.exists("/apps/ubuntu/24.04/bin/chrome")) {
  Sys.setenv(CHROMOTE_CHROME = "/apps/ubuntu/24.04/bin/chrome")
}

library(tidyverse)
library(rgl)

# ---- Distance functions (from .VAE_evo_model_vis.R) ----

calc_dist_along <- function(zseqs) {
  zdiffs <- diff(as.matrix(zseqs))
  sum(sqrt(rowSums(zdiffs^2)))
}

calc_rate_along <- function(zseqs, time_seqs) {
  z_dist <- calc_dist_along(zseqs)
  z_dist / (time_seqs[length(time_seqs)] - time_seqs[1])
}

calc_straight_dist <- function(zseqs) {
  m <- as.matrix(zseqs)
  sqrt(sum((m[nrow(m), ] - m[1, ])^2))
}

# ---- Color blending for trophic transitions ----

blend_colors <- function(col1, col2, t) {
  rgb1 <- col2rgb(col1) / 255
  rgb2 <- col2rgb(col2) / 255
  rgb((1 - t) * rgb1[1] + t * rgb2[1],
      (1 - t) * rgb1[2] + t * rgb2[2],
      (1 - t) * rgb1[3] + t * rgb2[3])
}

#' Compute smoothed colors for a trophic niche sequence
#' Blends over a window around transition points
compute_smooth_colors <- function(prtroph_seq, niche_pal, blend_window = 5) {
  n <- length(prtroph_seq)
  colors <- niche_pal[prtroph_seq]

  # Find transition points
  transitions <- which(prtroph_seq[-1] != prtroph_seq[-n])

  for (tr in transitions) {
    col_before <- niche_pal[prtroph_seq[tr]]
    col_after <- niche_pal[prtroph_seq[tr + 1]]
    half_w <- blend_window %/% 2

    blend_start <- max(1, tr - half_w + 1)
    blend_end <- min(n, tr + half_w)

    for (j in blend_start:blend_end) {
      t <- (j - blend_start) / (blend_end - blend_start)
      colors[j] <- blend_colors(col_before, col_after, t)
    }
  }
  colors
}

# ---- Mesh generation (wraps existing pattern from .VAE_evo_model_vis.R:91) ----

#' Generate mesh from decoded beak code (64-dim)
#' @param prcode 1-row matrix of 64-dim decoded beak codes
#' @param scaling list with $means and $sds (from code_dat_means_sds.rds)
#' @param beak_mod loaded SDF model on CUDA
#' @param resolution marching cubes resolution (default 100 for animation)
make_mesh <- function(prcode, scaling, beak_mod, resolution = 100) {
  code_scaled <- (prcode * scaling$sds) + scaling$means
  beak_mod$get_mesh(
    torch::torch_tensor(code_scaled, device = "cuda"),
    resolution = resolution,
    smooth = FALSE
  )
}

# ---- Frame rendering ----

# ---- Standard beak rotations (from get_latent_samples.R) ----
# These orient the beak mesh into recognizable viewing angles

#' Apply a standard beak rotation
#' @param mesh rgl mesh3d object
#' @param view one of "side", "side2", "threequarter",
#'   "threequarter2", "top", "front"
rotate_beak <- function(mesh, view = "side") {
  switch(view,
    side = mesh |>
      rotate3d(pi / 2, 0, 0, 1) |>
      rotate3d(pi / 2, 0, 1, 0),
    side2 = mesh |>
      rotate3d(-pi / 2, 0, 0, 1) |>
      rotate3d(-pi / 2, 0, 1, 0),
    threequarter = mesh |>
      rotate3d(pi / 2, 0, 0, 1) |>
      rotate3d(pi / 2, 0, 1, 0) |>
      rotate3d(3 * pi / 4, 0, 0, 1),
    threequarter2 = mesh |>
      rotate3d(pi / 2, 0, 0, 1) |>
      rotate3d(pi / 2, 0, 1, 0) |>
      rotate3d(pi / 4, 0, 0, 1),
    top = mesh |>
      rotate3d(pi, 1, 0, 0),
    front = mesh |>
      rotate3d(pi / 2, 0, 0, 1) |>
      rotate3d(pi / 2, 0, 1, 1),
    # default: side view
    mesh |>
      rotate3d(pi / 2, 0, 0, 1) |>
      rotate3d(pi / 2, 0, 1, 0)
  )
}

#' Render a single mesh to a PNG file
#' Uses webshot2 + Chrome for headless rendering
#' @param mesh rgl mesh3d object
#' @param col color string
#' @param filename output PNG path
#' @param width,height image dimensions
#' @param view beak viewing angle (see rotate_beak)
render_mesh_frame <- function(mesh, col, filename,
                               width = 800, height = 600,
                               view = "side",
                               theta = 0, phi = 0, fov = 30,
                               transparent_bg = FALSE) {
  # For "3/4" view, use raw mesh + camera angles
  # matching the published paper's perspective
  if (view == "3/4") {
    theta <- -150
    phi <- 20
    view <- "raw"
  }

  open3d()
  par3d(windowRect = c(0, 0, width, height))
  bg3d("white")

  # Add front-facing light to reduce dark shadows
  light3d(theta = -150, phi = 30, diffuse = "grey80",
          specular = "grey40")

  if (view == "raw") {
    shade3d(mesh, col = col)
  } else {
    rotated <- rotate_beak(mesh, view = view)
    shade3d(rotated, col = col)
  }

  view3d(theta = theta, phi = phi, fov = fov)
  # Use webshot2 (default) for headless rendering via Chrome
  snapshot3d(filename, width = width, height = height)
  close3d()

  # Make background transparent if requested
  if (transparent_bg) {
    img <- magick::image_read(filename)
    img <- magick::image_transparent(img, "white", fuzz = 5)
    magick::image_write(img, filename, format = "png")
  }
}

# ---- Timeline overlay ----

#' Add a timeline progress bar to a rendered frame
#' @param img magick image object
#' @param progress fraction 0-1 along the branch
#' @param width image width in pixels
#' @param y_pos y position from bottom (pixels)
#' @param bar_color color of the progress marker
#' @param label optional text to display near marker
add_timeline <- function(img, progress, width = 800,
                          y_pos = 60, bar_color = "grey30",
                          marker_color = "red",
                          label = NULL) {
  # Draw timeline bar and marker using magick
  margin <- round(width * 0.1)
  bar_width <- width - 2 * margin
  bar_y <- image_info(img)$height - y_pos

  # Create overlay with timeline
  overlay <- image_blank(
    image_info(img)$width,
    image_info(img)$height,
    color = "none"
  )

  # Draw the track line
  overlay <- image_draw(overlay)
  lines(
    c(margin, margin + bar_width),
    c(bar_y, bar_y),
    col = bar_color, lwd = 2
  )
  # Start and end ticks
  lines(c(margin, margin), c(bar_y - 6, bar_y + 6),
        col = bar_color, lwd = 2)
  lines(c(margin + bar_width, margin + bar_width),
        c(bar_y - 6, bar_y + 6),
        col = bar_color, lwd = 2)
  # Marker dot
  marker_x <- margin + round(bar_width * progress)
  points(marker_x, bar_y, pch = 19, cex = 2,
         col = marker_color)
  dev.off()

  # Composite overlay onto image
  img <- image_composite(img, overlay)

  # Add label if provided
  if (!is.null(label)) {
    img <- image_annotate(
      img, label, size = 14, color = bar_color,
      gravity = "south",
      location = paste0("+0+", y_pos + 10),
      font = "Helvetica"
    )
  }

  img
}

# ---- Arc-length re-parameterization ----

#' Resample a trajectory so frames are equally spaced
#' in Euclidean distance (arc length).
#' This prevents jerky animations caused by non-uniform
#' decoder expansion of latent space steps.
#' @param codes matrix (n_points x n_dims)
#' @param n_out number of output frames
#' @param troph character vector of trophic labels (same length as nrow(codes))
#' @return list with $codes (n_out x n_dims matrix) and $troph (character)
resample_by_arclength <- function(codes, n_out = nrow(codes),
                                   troph = NULL) {
  n <- nrow(codes)
  # Compute cumulative arc length
  diffs <- diff(codes)
  seg_lengths <- sqrt(rowSums(diffs^2))
  cum_length <- c(0, cumsum(seg_lengths))
  total_length <- cum_length[n]

  # Target equally-spaced arc lengths
  target <- seq(0, total_length, length.out = n_out)

  # Interpolate each dimension
  new_codes <- matrix(0, nrow = n_out, ncol = ncol(codes))
  for (d in seq_len(ncol(codes))) {
    new_codes[, d] <- approx(
      cum_length, codes[, d],
      xout = target, rule = 2
    )$y
  }

  # Interpolate trophic labels by nearest original frame
  new_troph <- NULL
  if (!is.null(troph)) {
    nearest_idx <- sapply(target, function(t) {
      which.min(abs(cum_length - t))
    })
    new_troph <- troph[nearest_idx]
  }

  list(codes = new_codes, troph = new_troph)
}

# ---- GIF assembly ----

#' Assemble PNG frames into animated GIF
#' @param png_files character vector of PNG file paths
#' @param gif_file output GIF path
#' @param fps frames per second
#' @param width,height GIF dimensions
assemble_gif <- function(png_files, gif_file, fps = 15,
                          width = 800, height = 600) {
  gifski::gifski(png_files, gif_file = gif_file,
                 delay = 1 / fps, width = width, height = height)
  message("GIF saved to: ", gif_file)
}

# ---- Trophic palette reconstruction ----

#' Reconstruct niche_pal without needing targets store
#' @param trophic_levels character vector of trophic niche level names
make_niche_pal <- function(trophic_levels) {
  pal <- Polychrome::createPalette(
    length(trophic_levels),
    wesanderson::wes_palettes$FantasticFox1[-1]
  )
  # createPalette gives generic names; replace with actual levels
  pal <- unname(pal)
  names(pal) <- trophic_levels
  pal
}

# ---- Try loading niche_pal from targets, fall back to reconstruction ----

load_niche_pal <- function() {
  pal <- tryCatch(
    targets::tar_read(niche_pal),
    error = function(e) NULL
  )

  # If loaded from targets, ensure names match trophic levels
  if (!is.null(pal)) {
    trophic_dat <- read_rds("data/trophic_dat.rds")
    trophic_levs <- levels(trophic_dat$trophic_niche)
    if (!all(trophic_levs %in% names(pal))) {
      names(pal) <- trophic_levs
    }
    return(pal)
  }

  # Fallback: reconstruct from scratch
  trophic_dat <- read_rds("data/trophic_dat.rds")
  trophic_levs <- levels(trophic_dat$trophic_niche)
  make_niche_pal(trophic_levs)
}
