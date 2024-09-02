library(reticulate)
library(tidyverse)

use_condaenv("ganspace")

torch <- import("torch")
np <- import("numpy")

latent_codes <- torch$load("data/model/sdf_net_latent_codes.to")$cpu()$detach()$numpy()
write_rds(latent_codes, "data/model/sdf_net_latent_codes.rds")

dat_pnts <- torch$load("data/model/sdf_points.to")$cpu()$numpy()
write_rds(dat_pnts, "data/model/sdf_points.rds")

dat_vals <- torch$load("data/model/sdf_values.to")$cpu()$numpy()
write_rds(dat_vals, "data/model/sdf_values.rds")

models <- list.files("data/model/checkpoints", 
                     full.names = TRUE)
latent_codes <- grep("sdf_net_latent_codes-epoch", models, value = TRUE)
models <- grep("sdf_net-epoch", models, value = TRUE)

weights <- map(models,
               ~ torch$load(.x) |>
                 map(function(x) x$cpu()$numpy()),
               .progress = TRUE)

files <- gsub("checkpoints", "checkpoints_rds", models)
files <- gsub(".to", ".rds", files)
walk2(weights, files,
      ~ write_rds(.x, .y),
      .progress = TRUE)

codes <- map(latent_codes,
             ~ torch$load(.x)$cpu()$detach()$numpy(),
             .progress = TRUE)

files <- gsub("checkpoints", "checkpoints_codes_rds", latent_codes)
files <- gsub(".to", ".rds", files)
walk2(codes, files,
      ~ write_rds(.x, .y),
      .progress = TRUE)

#mod1 <- torch$load("../data/model/checkpoints/sdf_net-epoch-00000.to")

#mod <- map(mod1, ~ .x$cpu()$numpy())

