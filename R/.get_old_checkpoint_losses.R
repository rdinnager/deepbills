library(fibre)
library(tidyverse)
library(torch)
library(cli)

sdfnet <- sdf_net()$cuda()$eval()
sdfnet$parameters %>% purrr::walk(function(param) param$requires_grad_(FALSE))

sdf_vals <- read_rds("data/model/sdf_values.rds")
sdf_points <- read_rds("data/model/sdf_points.rds")
latent_codes <- matrix(rnorm(2021 * 64, 0, 0.0001), nrow = 2021, ncol = 64)

chkpnts <- list.files("data/model/checkpoints_rds", full.names = TRUE)
chkpnts_codes <- list.files("data/model/checkpoints_codes_rds", full.names = TRUE)
chkpnts <- chkpnts[-1]
chkpnts_codes <- chkpnts_codes[-1]

nums <- c(1, 122, seq(150, 2450, by = 50), 41, 91, 20)
checkpoints <- chkpnts[nums]
checkpoint_codes <- chkpnts_codes[nums]

batch_size <- length(sdf_vals) / 2021
batch <- as.integer(as.numeric(gl(2021, batch_size)))

sdf_dataset <- dataset(name = "bb_ds",
                           initialize = function(pnts, vals, batch, latent_codes) {
                             self$pnts <- torch_tensor(pnts)
                             self$vals <- torch_tensor(vals) 
                             self$batch <- torch_tensor(batch)
                             self$latent_codes <- torch_tensor(latent_codes)
                           },
                           .getbatch = function(i) {
                             list(pnts = self$pnts[i, ], vals = self$vals[i], codes = self$latent_codes[self$batch[i], ])
                           },
                           replace_codes = function(codes) {
                             self$latent_codes <- torch_tensor(codes)
                           },
                           .length = function() {
                             self$pnts$size()[[1]]
                           })
train_ds <- sdf_dataset(sdf_points, sdf_vals, batch, latent_codes)
train_dl <- dataloader(train_ds, 300000, shuffle = TRUE)

rm(sdf_vals, sdf_points, latent_codes, batch)
gc()

#test <- train_dl$.iter()$.next()

# batchnum <- 0
# batchlen <- length(train_dl)
# recons <- numeric(batchlen)
# coro::loop(for (b in train_dl) {
#   batchnum <- batchnum + 1
#   test_out <- sdfnet(b$pnts$cuda(), b$codes$cuda())$clamp(-0.1, 0.1)
#   recon_loss <- mean(abs(b$vals$cuda()$clamp(-0.1, 0.1) - test_out))
#   recons[batchnum] <- as.numeric(recon_loss$cpu())
#   print(batchnum / batchlen)
# })
# 
# init_recon <- mean(recons)
# write_rds(init_recon, "data/model/losses/sdf_net-epoch-00000.rds")

load_weights <- function(sdfnet, sdf_weights) {
  with_no_grad({
    sdfnet$layers1$`0`$weight$set_(torch_tensor(sdf_weights$layers1.0.weight)$cuda())
    sdfnet$layers1$`0`$bias$set_(torch_tensor(sdf_weights$layers1.0.bias)$cuda())
    sdfnet$layers1$`2`$weight$set_(torch_tensor(sdf_weights$layers1.2.weight)$cuda())
    sdfnet$layers1$`2`$bias$set_(torch_tensor(sdf_weights$layers1.2.bias)$cuda())
    sdfnet$layers1$`4`$weight$set_(torch_tensor(sdf_weights$layers1.4.weight)$cuda())
    sdfnet$layers1$`4`$bias$set_(torch_tensor(sdf_weights$layers1.4.bias)$cuda())
    sdfnet$layers2$`0`$weight$set_(torch_tensor(sdf_weights$layers2.0.weight)$cuda())
    sdfnet$layers2$`0`$bias$set_(torch_tensor(sdf_weights$layers2.0.bias)$cuda())
    sdfnet$layers2$`2`$weight$set_(torch_tensor(sdf_weights$layers2.2.weight)$cuda())
    sdfnet$layers2$`2`$bias$set_(torch_tensor(sdf_weights$layers2.2.bias)$cuda())
    sdfnet$layers2$`4`$weight$set_(torch_tensor(sdf_weights$layers2.4.weight)$cuda())
    sdfnet$layers2$`4`$bias$set_(torch_tensor(sdf_weights$layers2.4.bias)$cuda())
  })
  
  sdfnet
}

mod_file <- checkpoints[1]
code_file <- checkpoint_codes[1]
recon_loss <- function(mod_file, code_file) {
  sdf_weights <- read_rds(mod_file)
  sdfnet <- load_weights(sdfnet, sdf_weights)
  latent_codes <- read_rds(code_file)
  train_ds$replace_codes(latent_codes)
  train_dl <- dataloader(train_ds, 500000, shuffle = TRUE)
  batchnum <- 0
  batchlen <- length(train_dl)
  recons <- numeric(batchlen)
  coro::loop(for (b in train_dl) {
    batchnum <- batchnum + 1
    test_out <- sdfnet(b$pnts$cuda(), b$codes$cuda())$clamp(-0.1, 0.1)
    recon_loss <- mean(abs(b$vals$cuda()$clamp(-0.1, 0.1) - test_out))
    recons[batchnum] <- as.numeric(recon_loss$cpu())
    cli_progress_message(batchnum / batchlen)
  })
  recon <- mean(recons)
  save_file <- gsub("checkpoints_rds", "losses", mod_file)
  write_rds(recon, save_file)
  recon
}

recon_losses <- map2(checkpoints, checkpoint_codes, recon_loss,
                     .progress = TRUE)

