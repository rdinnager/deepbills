library(reticulate)
library(readr)
library(dplyr)
library(ggplot2)
library(directlabels)
library(plotly)
library(purrr)
library(rgl)
library(ape)
library(mapdeck)
library(furrr)
library(RRphylo)

reticulate::use_condaenv("tf-gpu", required = TRUE)

np <- import("numpy")
torch <- import("torch")
umap <- import("umap")

meta <- readr::read_csv("data/resource.csv")
taxo <- readr::read_csv("data/Jetz/2012-03-04206D-master_taxonomy.csv")
#file_names <- readr::read_lines("Marian_results/filenames.txt")
file_names <- list.files("data/sdf")

tree <- ape::read.tree("data/phylogenies/Stage2_MayrParSho_Ericson_set1_decisive.tre")

latent_codes <- (torch$load("models/sdf_net_latent_codes.to") %>%
  torch$Tensor$cpu())$detach()$numpy()

meta_taxo <- meta %>%
  dplyr::left_join(taxo,
                   by = c("Binomal_Jetz" = "TipLabel"))

fams <- dplyr::tibble(ID = gsub(".obj", "", file_names)) %>%
  dplyr::left_join(meta_taxo %>%
                     dplyr::select(ID, family = BLFamilyEnglish)) %>%
  dplyr::mutate(fam_num = as.numeric(as.factor(family)))

umapper <- umap$UMAP(min_dist = 0.75)
umapper3d <- umap$UMAP(min_dist = 0.75, n_components = 3L)

umap_latent <- umapper$fit_transform(latent_codes, fams$fam_num)
plot(umap_latent)

umap_latent3d <- umapper3d$fit_transform(latent_codes, fams$fam_num)

dat <- meta_taxo %>%
  dplyr::right_join(dplyr::tibble(ID = gsub(".obj", "", file_names),
                                  X = umap_latent[ , 1],
                                  Y = umap_latent[ , 2],
                                  X0 = umap_latent3d[ , 1],
                                  Y0 = umap_latent3d[ , 2],
                                  Z0 = umap_latent3d[ , 3]))

unique(dat$BLFamilyEnglish)

# pdf("figures/test_umap.pdf", height = 30, width = 30)
# 
# ggplot(dat, aes(X, Y)) +
#   geom_point(aes(colour = BLFamilyEnglish), alpha = 0.7, size = 6) +
#   geom_dl(aes(label = BLFamilyEnglish, colour = BLFamilyEnglish), method = "ahull.grid") +
#   scale_color_discrete(guide = "none") +
#   theme_minimal()
# 
# dev.off()

p <- ggplot(dat, aes(X, Y)) +
  geom_point(aes(colour = BLFamilyEnglish), alpha = 0.5) +
  scale_color_discrete(guide = "none") +
  theme_minimal()

pl <- plotly::ggplotly(p)
pl

dat <- dat %>%
  dplyr::mutate(latent_code = purrr::map(1:nrow(latent_codes), ~ latent_codes[.x, ]))


umap_latent2 <- umapper$fit_transform(latent_codes)
plot(umap_latent2)

dat2 <- meta_taxo %>%
  dplyr::right_join(dplyr::tibble(ID = gsub(".obj", "", file_names),
                                  X = umap_latent2[ , 1],
                                  Y = umap_latent2[ , 2]))

dat2 <- dat2 %>%
  dplyr::mutate(latent_code = purrr::map(1:nrow(latent_codes), ~ latent_codes[.x, ]))

p <- ggplot(dat2, aes(X, Y)) +
  geom_point(aes(colour = BLFamilyEnglish), alpha = 0.5) +
  scale_color_discrete(guide = "none") +
  theme_minimal()

pl <- plotly::ggplotly(p)
pl


source("R/sdf_tools.R")
setup_SDF()

test <- get_meshes_from_latent(latent_codes[100, ], show = TRUE)
test <- get_meshes_from_latent(latent_codes[1845, ], show = TRUE, voxel_res = 256L)

test <- get_meshes_from_latent(dat$latent_code[[4]], show = TRUE)
test2 <- get_meshes_from_latent(dat$latent_code[[100]], show = TRUE, voxel_res = 128L)
test3 <- get_meshes_from_latent(dat$latent_code[[800]], show = TRUE, voxel_res = 128L)

dat <- dat %>%
  dplyr::mutate(recon_mesh = pbapply::pblapply(latent_code, get_meshes_from_latent, voxel_res = 128L))

rgl::shade3d(dat$recon_mesh[[100]], col = "#d2b232")

readr::write_rds(dat, "data/latent_code_reconstructions.rds", compress = "gz")

future::plan(future::multiprocess())

# dat <- dat %>%
#   dplyr::mutate(recon_mesh_256 = furrr::future_map(latent_code, 
#                                                    ~get_meshes_from_latent(.x, voxel_res = 256L),
#                                                    .progress = TRUE))

dat <- dat %>%
  dplyr::mutate(recon_mesh_256 = pbapply::pblapply(latent_code, get_meshes_from_latent, voxel_res = 256L))
readr::write_rds(dat, "data/latent_code_reconstructions2.rds", compress = "gz")

rgl::shade3d(dat$recon_mesh_256[[100]], col = "#d2b232")

dat <- readr::read_rds("data/latent_code_reconstructions.rds")

test <- get_meshes_from_latent(dat$latent_code[[10]], voxel_res = 256L)
rgl::shade3d(test, col = "#d2b232")
