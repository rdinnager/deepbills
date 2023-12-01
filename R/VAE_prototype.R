####### Load Libraries ###################
library(targets)
library(torch)
library(tidyverse)
library(dagnn) ## This is in-development package: https://github.com/rdinnager/dagnn
library(zeallot)
library(conflicted)
library(phyf)
library(torchopt)

conflicts_prefer(torch::optim_adamw)

####### Load Data ########################

avonet <- tar_read(bird_beak_avonet)

code_dat <- avonet |>
  select(starts_with("latent_code_")) |>
  as.matrix() |>
  scale()

write_rds(list(means = attr(code_dat, "scaled:center"),
               sds = attr(code_dat, "scaled:scale")),
          "data/code_dat_means_sds.rds")

trophic_dat <- avonet |>
  select(label, Trophic.Niche) |>
  mutate(trophic_niche = as.factor(Trophic.Niche))

write_rds(trophic_dat, "data/trophic_dat.rds")

trophic_one_hot <- trophic_dat |>
  select(label, trophic_niche) |>
  mutate(yes = 1) |>
  pivot_wider(names_from = trophic_niche, values_from = yes, values_fill = 0)

trophic_one_hot <- trophic_one_hot |>
  select(-label) |>
  as.matrix()

trophic_dat <- trophic_dat |>
  select(trophic_niche) |>
  pull(trophic_niche)

avonet_dataset <- dataset(name = "avonet_ds",
                          initialize = function(codes, trophic_in, trophic_out) {
                            self$codes <- torch_tensor(codes)
                            self$trophic_in <- torch_tensor(trophic_in)
                            self$trophic_out <- torch_tensor(trophic_out)
                          },
                          .getbatch = function(i) {
                            list(codes = self$codes[i, ],
                                 trophic_in = self$trophic_in[i],
                                 trophic_out = self$trophic_out[i, ])
                          },
                          .length = function() {
                            self$codes$size()[[1]]
                          })

train_ds <- avonet_dataset(code_dat, trophic_dat, trophic_one_hot)
train_dl <- dataloader(train_ds, 506, shuffle = TRUE)

b <- train_dl$.iter()$.next()

####### Setup VAE model as nn_module ###########
vae_mod <- nn_module("CVAE",
                 initialize = function(codes_dim, trophic_dim, trophic_embed_dim, latent_dim, breadth) {
                   self$latent_dim <- latent_dim
                   self$codes_dim <- codes_dim
                   self$trophic_dim <- trophic_dim
                   self$trophic_embed_dim <- trophic_embed_dim
                   self$trophic_embedding <- nn_embedding(trophic_dim, trophic_embed_dim)
                   input_dim <- codes_dim + trophic_embed_dim
                   self$input_dim <- input_dim
                   self$encoder <- nndag(x_1 = ~ input_dim,
                                         e_1 = x_1 ~ breadth,
                                         e_2 = e_1 ~ breadth,
                                         e_3 = e_2 ~ breadth,
                                         means = e_3 ~ latent_dim,
                                         logvars = e_3 ~ latent_dim,
                                         .act = list(nn_relu,
                                                     logvars = nn_identity,
                                                     means = nn_identity))

                   self$decoder <- nndag(z_1 = ~ latent_dim,
                                         d_1 = z_1 ~ breadth,
                                         d_2 = d_1 ~ breadth,
                                         d_3 = d_2 ~ breadth,
                                         out_codes = d_3 ~ codes_dim,
                                         out_trophic = d_3 ~ trophic_dim,
                                         .act = list(nn_relu,
                                                     out_codes = nn_identity,
                                                     out_trophic = nn_identity))

                 },
                 reparameterize = function(mean, logvar) {
                   std <- torch_exp(torch_tensor(0.5, device = "cuda") * logvar)
                   eps <- torch_randn_like(std)
                   eps * std + mean
                 },
                 loss_function = function(recon_codes, input_codes, recon_trophic, input_trophic, mean, log_var,  w_tr, loggamma) {
                   kl <- torch_sum(torch_exp(log_var) + torch_square(mean) - log_var, dim = 2L) - self$latent_dim
                   recon1 <- torch_sum(torch_square(input_codes - recon_codes), dim = 2L) / torch_exp(loggamma)
                   recon2 <- nnf_cross_entropy(recon_trophic, input_trophic, reduction = "none") / torch_exp(loggamma)
                   loss <- torch_mean(recon1 / 2 + w_tr * recon2 / 2 + kl)
                   list(loss, torch_mean(recon1) * torch_exp(loggamma), torch_mean(recon2) * torch_exp(loggamma), torch_mean(kl))
                 },
                 forward = function(codes, trophic_in) {
                   tr <- self$trophic_embedding(trophic_in)
                   c(means, log_vars) %<-% self$encoder(torch_cat(list(codes, tr), dim = 2L))
                   z <- self$reparameterize(means, log_vars)
                   recon <- self$decoder(z)
                   
                   list(codes, recon, means, log_vars)
                 })

######## Run VAE stage 1 training ##################
## params
codes_dim <- ncol(code_dat)
trophic_dim <- ncol(trophic_one_hot)
trophic_embed_dim <- 16L
latent_dim <- 64L
breadth <- 1024L

vae <- vae_mod(codes_dim, trophic_dim, trophic_embed_dim, latent_dim, breadth)
vae <- vae$cuda()

num_epochs <- 50000

lr <- 0.002
optimizer <- optim_adam(vae$parameters, lr = lr)
scheduler <- lr_one_cycle(optimizer, max_lr = lr,
                          epochs = num_epochs, steps_per_epoch = length(train_dl),
                          cycle_momentum = FALSE)

gamma_x <- 0.002
mseloss <- 0.002#gamma_x^2
loggamma <- log(gamma_x)

w_tr <- 1 / 64

for (epoch in 1:num_epochs) {

    batchnum <- 0
    coro::loop(for (b in train_dl) {

        batchnum <- batchnum + 1
        optimizer$zero_grad()
        c(input_codes, reconstruction, mean, log_var) %<-% vae(b$codes$cuda(), b$trophic_in$cuda())
        c(loss, recon_code_loss, recon_trophic_loss, kl_loss) %<-%
          vae$loss_function(reconstruction$out_codes, input_codes, 
                            reconstruction$out_trophic, b$trophic_in$cuda(), 
                            mean, log_var, w_tr,
                            torch_tensor(loggamma, device = "cuda")
                            )

        mseloss <- min(mseloss, mseloss * .99 + as.numeric(((recon_code_loss$cpu() / 2) + (recon_trophic_loss$cpu() / 2)) / vae$input_dim) * .01)
        gamma_x <- sqrt(mseloss)
        loggamma <- log(gamma_x)

        if(batchnum %% 2 == 0) {

            cat("Epoch: ", epoch,
                "    batch: ", batchnum,
                "    loss: ", as.numeric(loss$cpu()),
                "    code recon loss: ", as.numeric(recon_code_loss$cpu()),
                "    trophic recon loss: ", as.numeric(recon_trophic_loss$cpu()),
                "    KL loss: ", as.numeric(kl_loss$cpu()),
                "    loggamma: ", as.numeric(loggamma),
                #"    loggamma: ", loggamma,
                "    active dims: ", as.numeric((torch_exp(log_var)$mean(dim = 1L) < 0.5)$sum()$cpu()),
                "\n")

        }
        loss$backward()

        # if((epoch / num_epochs) > 0.6) {
        #   nn_utils_clip_grad_norm_(vae$parameters, 1)
        # }

        optimizer$step()
        scheduler$step()
    })
}

options(torch.serialization_version = 2)
torch_save(vae, "data/bill_vae_w_trophic_v1.to")

vae <- torch_load("data/bill_vae_w_trophic_v1.to")
vae <- vae$cuda()

post_dl <- dataloader(train_ds, 506, shuffle = FALSE)
post_it <- as_iterator(post_dl)

all_mean_var <- list()
it <- 0
loop(for(i in post_dl) {
  it <- it + 1
  with_no_grad({
    all_mean_var[[it]] <- vae(i$codes$cuda(), i$trophic_in$cuda())
  })
  print(it)
})

all_means <- map(all_mean_var, 3) %>%
  map(~ as.matrix(.x$cpu())) %>%
  do.call(rbind, .)

all_vars <- map(all_mean_var, 4) %>%
  map(~ as.matrix(torch_exp(.x)$cpu())) %>%
  do.call(rbind, .)

vars <- apply(all_vars, 2, mean)

vars_df <- tibble(`Mean Variance` = vars) %>%
  mutate(Status = ifelse(`Mean Variance` < 0.5, "Active", "Ignored"))

ggplot(vars_df, aes(`Mean Variance`)) +
  geom_histogram(aes(fill = Status)) +
  ylab("Count") +
  theme_minimal() +
  theme(legend.position = c(0.5, 0.5))

active_dims <- which(vars < 0.5)
plot(all_means[ , active_dims][ , 1:2])
rgl::points3d(all_means)

write_rds(active_dims, "data/active_dims_16dim.rds")

######### Save latent codes ################

active_dims <- which(vars < 0.5)
latent_means <- all_means[ , active_dims]
latent_vars <- all_vars[ , active_dims]

colnames(latent_means) <- colnames(latent_vars) <- paste0("latent_", 1:ncol(latent_means))

latent_df <- as.data.frame(latent_means) %>%
  mutate(Species = avonet$label) %>%
  pivot_longer(-Species, names_to = "latent_dim", values_to = "latent_mean") %>%
  left_join(as.data.frame(latent_vars) %>%
              mutate(Species = avonet$label) %>%
              pivot_longer(-Species, names_to = "latent_dim", values_to = "latent_var")) %>%
  mutate(latent_lower = latent_mean - 1.96 * sqrt(latent_var),
         latent_upper = latent_mean + 1.96 * sqrt(latent_var))

write_csv(latent_df, "data/bills_vae_latent_codes_16dim.csv")
latent_df_means <- as.data.frame(latent_means) %>%
  mutate(Species = avonet$label) %>%
  select(Species, everything())

write_csv(latent_df_means, "data/bills_vae_latent_codes_means_only_16dim.csv")

