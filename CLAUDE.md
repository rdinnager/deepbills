# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Deep Bills is a research project analyzing 3D bird bill (beak) morphology using deep learning and phylogenetic comparative methods. The published work ([Dinnage & Kleineberg 2025, PLoS Computational Biology](https://doi.org/10.1371/journal.pcbi.1012887)) trained a DeepSDF model on 3D scans of bills from 2,020 bird species (from the NHM Mark My Bird dataset), learning 64-dimensional latent representations that predict trophic niche with high accuracy (balanced accuracy ~0.79) and capture stronger phylogenetic signal than traditional PCA on landmarks (Blomberg's K = 0.81 vs 0.57). A conditional VAE generates realistic beak shapes conditioned on trophic niche.

The current active development direction models beak shape evolution on the **Riemannian manifold** implied by the VAE latent space, using non-Euclidean distances rather than assuming the latent space is flat. This is the `.VAE_*` script series described below. The latest iteration is the **v3 Bayesian model** (`.VAE_evo_model_v3_bayesian.R`), which switches the manifold distance from path **length** (arc length, `sqrt(vᵀGv)`) to path **energy** (the squared Riemannian norm `vᵀGv`, no sqrt). This makes the whole fit a proper **Bayesian MAP estimate of Brownian motion on the manifold** — the Onsager–Machlup action for true BM — with every loss weight an interpretable precision. The full derivation is in [`notes/bayesian_riemannian_evolution.md`](notes/bayesian_riemannian_evolution.md).

**Open scientific question (2026-07-10, from a meeting with Shinichi Nakagawa):** the model can reconstruct full beak shapes on the phylogeny, but the paper needs a compelling biological hypothesis to anchor it, not just a capable tool. Getting Shinichi's/Azumi's input on that framing is the current priority. See [`../_notes/meetings/2026-07-10-shinichi.md`](../_notes/meetings/2026-07-10-shinichi.md).

## Build & Run Commands

This is a `targets`-based R pipeline project (not a package). There is no traditional build/test/lint system.

```r
# Run the full pipeline (skips up-to-date targets)
targets::tar_make()

# Visualize the dependency graph
targets::tar_visnetwork()

# Check which targets are outdated
targets::tar_outdated()

# Load a specific target result into the R session
targets::tar_read(target_name)

# Run a specific target and its upstream dependencies
targets::tar_make(names = target_name)

# Invalidate and rebuild a target
targets::tar_invalidate(target_name)
```

All R source files in `R/` are auto-sourced at pipeline start via `lapply(list.files("./R", full.names = TRUE), source)` in `_targets.R`.

## Architecture

### Pipeline (`_targets.R`)

The ~1200-line `_targets.R` defines the entire analysis as a DAG of targets. Major stages:

1. **Data ingestion** - AVONET trait data, phylogenetic trees (ClaDS), landmark coordinates, diet databases
2. **Deep learning** - Train VAE and CVAE models on beak morphology (`run_vae`, `run_cvae`), encoding 3D shape into 64-dimensional latent codes
3. **Latent extraction** - Extract latent codes from trained models (`get_vae_latents`, `get_cvae_latents`)
4. **Gaussianization** - Lambert W transform to normalize heavy-tailed latent distributions (`gaussianize`)
5. **Phylogenetic modeling** - Brownian motion evolutionary models on latent codes using `fibre` + INLA
6. **Trophic analysis** - Random forest classification of trophic niches from latent codes, cross-validation
7. **Visualization** - Latent space sampling, trajectory images, UMAP plots, phylogenetic tree figures

### Key Source Files (`R/`)

- `run_cvae.R` / `run_vae.R` - Train conditional/unconditional VAE models using `torch` and `dagnn`
- `get_cvae_latents.R` / `get_vae_latents.R` - Extract latent representations from trained models
- `get_latent_samples.R` - Sample latent space and render 3D beak reconstructions
- `gaussianize.R` - Lambert W transformation for normalizing latent codes
- `reconcile_clads.R` - Join ClaDS diversification rates with morphological data
- `fit_RF_trophic.R` / `tune_RF_trophic.R` - Random forest trophic niche classification
- `make_diet_data.R` (called via `make_diet_data()`) - Process avian diet database
- `landmarks_functions.R` - Landmark-based morphometric analysis utilities
- `functions.R` - General helpers (e.g., `get_segments` for torch tensor operations)
- Files prefixed with `.` (e.g., `.VAE_evo_model.R`) are standalone scripts for the Riemannian manifold evolutionary modeling (see below)

### Python Subproject (`shapegan/`)

Adapted from [ShapeGAN](https://github.com/ndoll/ShapeGAN) for 3D shape generation. Contains training scripts for autoencoders, GANs (WGAN), and SDF-based models on voxelized beak data. `prepare_data_birds.py` handles bird-specific data prep.

### Data (`data/` - gitignored)

Large datasets stored locally, not in version control. Key files:
- `AVONET3_BirdTree.csv` - Avian morphological traits
- `AvianDietDatabase.txt` - Diet/trophic data
- `landmarks_oriented/` - 3D beak landmark coordinates
- `clads_pf.rds` - ClaDS phylogenetic diversification rates

## Key Technical Details

- **Torch serialization**: `options(torch.serialization_version = 2)` is set globally; torch targets use `format = "torch"`
- **Namespace conflicts**: Managed explicitly via `conflicted` package (`dplyr::filter`, `dplyr::select`, `imager::save.image`, `zeallot::%<-%`)
- **Pipeline error handling**: `tar_option_set(error = "continue")` - pipeline continues past individual target failures
- **Parallel execution**: Uses `future` package; `future.globals.onReference = "error"` catches non-serializable references
- **INLA models**: Run with `inla.mode = "experimental"` and typically 6 cores
- **Phylogenetic data**: Uses `phyf` phylogenetic frames (tibble-like structures with phylogenetic metadata)

## Riemannian Manifold Evolutionary Model (Active Development)

This is the current focus of the project: modeling beak shape evolution on the Riemannian manifold induced by the VAE, rather than treating the latent space as Euclidean. The `.VAE_*` scripts in `R/` implement this pipeline as standalone scripts (not yet integrated into `_targets.R`).

### Core Idea

A VAE's latent space has a natural Riemannian geometry: the posterior distribution defines a metric tensor that makes movement "cheap" near data-dense regions and "expensive" in unsupported regions. Evolutionary change along phylogenetic branches is modeled as curved paths through this manifold, with path lengths computed using the local metric rather than Euclidean distance.

### Pipeline Stages

```
.VAE_prototype.R          Train CVAE (64-dim latent, 16 active dims)
        |
.VAE_aces_init.R          Brownian motion ancestral states via fibre/INLA (initialization)
        |
.VAE_evo_model.R (v1)     Optimize rates + curved paths on Riemannian manifold
.VAE_evo_model_v2.0.R     v2: adds decoder-space losses (beak shape + trophic smoothness)
        |
.VAE_evo_post.R           Extract predictions (no_grad forward pass)
        |
.VAE_evo_prediction_analysis.R      Decode latent trajectories -> beak shapes + trophic niches
.VAE_evo_prediction_analysis_linear.R   Same but for straight-line baseline
        |
.VAE_evo_model_vis.R      Rate computation, 3D beak meshes, colored phylogenies
.VAE_evo_model_vis_v2.R   Taxonomic annotation, trophic transition matrices, polar trees
.VAE_tree_mesh_process.R   UMAP embedding of all latent trajectories
        |
.VAE_compare_trophic_model.R   Comparison with standard Mk discrete ancestral reconstruction
.phylo_pca_fibre.R              Comparison with direct landmark Brownian motion model
```

### Key Scripts in Detail

**`.VAE_prototype.R`** - Trains the conditional VAE using `torch`/`dagnn`. Architecture: encoder (3 x 1024 hidden layers) maps beak shape codes + one-hot trophic niche to 64-dim latent; decoder reconstructs both. Loss = KL divergence + MSE reconstruction + cross-entropy on trophic niche. 50,000 epochs, Adam + one-cycle LR. Post-training identifies 16 "active" dimensions (posterior variance < 0.5). Saves model as `bill_vae_w_trophic_v1.to`.

**`.VAE_aces_init.R`** - Computes initial ancestral state estimates via `fibre(latent_1 + ... + latent_16 ~ bre_brownian(phlo))` using INLA. Extracts per-edge evolutionary rates as starting values for the manifold optimizer. Saves `init_rates_16dim.csv`.

**`.VAE_evo_model.R`** (v1) - The core Riemannian optimization. Key components:

- **Metric tensor** (`get_metric_tensor`): `G(z)_jj = 1 / [sum_i (1/sigma_ij^2) * exp(-d_Mah(z, mu_i)^2 / rho^2) + lambda]`. Diagonal metric that is small near VAE centroids (cheap movement) and large far from data (expensive movement). `rho` = bandwidth, `lambda` = regularization floor.
- **Curved paths** (`get_segments`): Cubic polynomial parameterization `y(t) = a*t^3 + b*t^2 + (len-a-b)*t + z_start` with learnable `a`, `b` curvature parameters per edge/dimension.
- **Rho annealing**: Cosine schedule from `3 * max_min_dist` to `max_min_dist / 3`. Starts nearly Euclidean, progressively sharpens manifold structure (simulated annealing analog).
- **Loss**: manifold path length + tip_weight * tip MSE.
- **`mani_evo_mod`** (nn_module): Learnable parameters are `rates`, `a`, `b` per edge per dimension. Tree structure encoded via sparse phylogenetic matrices from `phyf`.

**`.VAE_evo_model_v2.0.R`** - Adds decoder-aware losses: runs latent path points through the frozen VAE decoder during training. Total loss = manifold distance + (1/64) * decoded beak shape change + (1/10) * decoded trophic niche change + 10 * tip loss + (1/100) * root loss. Adds a learnable `root_values` parameter.

**`.VAE_evo_model_v3_bayesian.R`** (latest, added on the `evoV3` branch 2026-07-10) - Reformulates the whole objective as a Bayesian MAP estimate. **Key change:** uses path **energy** `get_manifold_energy(vel, metric) = (vel*metric*vel)$sum(dim=2)` (squared Riemannian norm, no sqrt) instead of the old arc-length `get_manifold_dist` (kept for reference). The squared/energy form is the Onsager–Machlup action for Brownian motion on the manifold, so it corresponds to a true BM (Gaussian) prior, is smooth at v=0 (no gradient kink), and each loss weight maps to a precision (1/2σ²). Writes `data/v3_bayesian/mani_evo_mod_v3_bayesian.to` + checkpoints, and the ancestral estimates `bill_vae_aces_16dim_v3_bayesian(.linear).rds`. Companion `.VAE_evo_model_v3_extract.R` extracts the ancestral estimates from a trained model; `.animation_*.R` render evolutionary trajectories (divergence, edge-scan, UMAP phenogram). Full theory in `notes/bayesian_riemannian_evolution.md`.

**`.VAE_evo_prediction_analysis.R`** - Decodes optimized latent trajectories (50 time points per edge) through the VAE decoder to get predicted beak shapes and trophic niches at every point along every branch. `decode_zseqs()` embeds 16 active dims back to 64, runs decoder, extracts argmax trophic niche.

**`.VAE_compare_trophic_model.R`** - Baseline comparison using `castor::asr_mk_model()` (ARD rate matrix) for discrete ancestral trophic niche reconstruction.

**`.phylo_pca_fibre.R`** - Alternative comparison: Brownian motion model directly on PCA-whitened landmark coordinates (pre-whitened via phylogenetic precision matrix from `MCMCglmm::inverseA()`).

### Key Saved Artifacts

- `bill_vae_w_trophic_v1.to` - Trained CVAE model
- `active_dims_16dim.rds` - Which 16 of 64 latent dims are active
- `bills_vae_latent_codes_16dim.csv` - Per-species latent means + variances
- `init_rates_16dim.csv` - Initial Brownian motion rates per edge
- `mani_evo_mod_run_new_full_param_16dim_rho_scedule.to` - Trained v1 manifold model
- `bill_vae_aces_16dim_noise_schedule_v2.rds` - v2 curved predictions
- `bill_vae_aces_16dim_noise_schedule_v2_linear.rds` - v2 linear baseline predictions
- `evo_predictions_all_edges_v2.rds` - Decoded predictions for all edges

### Additional Dependencies for Manifold Pipeline

Beyond the base dependencies: `FNN` (k-NN for rho computation), `GPUmatrix` (sparse matrices to GPU torch tensors), `Matrix` (sparse tree structure), `castor` (Mk model comparison), `MCMCglmm` (phylogenetic precision matrix), `diagram` (transition matrix visualization).

## Dependencies

R packages are loaded in `packages.R`. Core dependencies: `targets`, `torch`, `dagnn`, `fibre`, `phyf`, `INLA`, `ape`, `tidyverse`, `tidymodels`, `Morpho`, `rgl`, `ggtree`, `uwot`.

## Git & hub integration (git-on-Drive scheme)

This project lives in Russell's hub (`G:\Shared drives\COBL Data\Projects\deepbills`).
The work tree is on Google Drive but the **`.git` lives off Drive** at
`C:\Users\dinnage\git-dirs\deepbills` (Drive corrupts a normal `.git`). Use the
`git-drive` wrapper, never plain `git` from the work tree:

```bash
cd "/g/Shared drives/COBL Data/Projects/_scratch/git-drive"
./git-drive deepbills <git args…>     # status, log, add, commit, push …
./git-drive sync deepbills            # fetch (the local mirror is machine-local)
```

Remote: `github.com/rdinnager/deepbills` (public). **Always `sync` before trusting
branch/commit state.** `data/`, `output/`, `_targets/`, `figures/`, and the `doc/`
manuscript materials are git-ignored and kept **local-only** (they never go to the
public repo); the bulk data + trained models live on Drive, backed up from the
`E:\Projects\deepbills` external-HDD copy.

## Progress

<!-- PROGRESS:BEGIN — maintained by the weekly /project-status scan; hand-edits welcome and respected -->
**Stage:** active (reviving toward publication)
**Last reviewed:** 2026-07-10
**Current status:** Assembled into the hub 2026-07-10: cloned from GitHub, backfilled
with the full 97 GB authoritative copy from `E:\Projects\deepbills`, and reconciled —
`main` now contains E:'s recovered 2024-25 work **and** the merged `evoV3` branch (the
v3 Bayesian path-energy model, extractor, animation scripts, and the 772-line theory
notes). All data + the trained v3 model + HPC outputs (`data/v3_bayesian/`, edge-scan,
animations) are present locally. The published Dinnage & Kleineberg 2025 (PLoS Comp
Biol) paper is the methods foundation; this repo is its evolutionary extension.

**Known / resolved — do NOT re-flag:**
- Missing data — RESOLVED 2026-07-10 (backfilled from the E: external HDD; kept local,
  git-ignored). Original heavy data lived on a Linux box (`/data/dinnage/Projects/`).
- The squared-distance evo model exists: it's the v3 Bayesian model on the merged
  `evoV3` work (path **energy** vs the old arc **length**).

**Open items:**
- **Generate the v3 ancestral estimates** — `bill_vae_aces_16dim_v3_bayesian.rds` and
  `..._linear.rds` were never created (the extract step wasn't run). Run
  `.VAE_evo_model_v3_extract.R` (all inputs are present locally; needs torch/GPU).
- **Write-up + biological question** — assemble a project description and get
  Shinichi's/Azumi's input on the anchoring hypothesis (see the meeting note).
- **Verify the pipeline runs** end-to-end (torch/CUDA; the `_targets` cache should let
  `tar_make()` skip most rebuilds; watch for hardcoded `/data/dinnage/...` Linux paths).
- Line-ending (CRLF) churn on many tracked files from the E: copy is left as-is;
  discard/normalize when convenient (content is identical, safe on the E: copy).
<!-- PROGRESS:END -->

