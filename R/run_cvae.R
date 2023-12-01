#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param bird_beak_avonet
run_cvae <- function(bird_beak_avonet) {

  codes <- bird_beak_avonet %>%
    select(starts_with("latent_")) %>%
    as.matrix() %>%
    scale()
  
  trophic <- bird_beak_avonet %>%
    select(Trophic.Niche) %>%
    mutate(Trophic.Niche = factor(Trophic.Niche)) %>%
    pull(Trophic.Niche)
  
  bb_dataset <- dataset(name = "bb_ds",
                           initialize = function(codes, trophic) {
                             self$codes <- torch_tensor(codes)
                             self$troph <- torch_tensor(trophic) 
                           },
                           .getbatch = function(i) {
                             list(codes = self$codes[i, ], troph = self$troph[i])
                           },
                           .length = function() {
                             self$codes$size()[[1]]
                           })
  train_ds <- bb_dataset(codes, trophic)
  train_dl <- dataloader(train_ds, 1012, shuffle = TRUE)

  vae_mod <- nn_module("CVAE",
                 initialize = function(input_dim, c_dim, c_embed_dim, latent_dim, breadth = 1024L, loggamma_init = 0) {
                   self$latent_dim <- latent_dim
                   self$input_dim <- input_dim
                   self$c_dim <- c_dim
                   self$c_embed_dim <- c_embed_dim
                   self$encoder <- nndag(y = ~ input_dim,
                                         c = ~ c_embed_dim,
                                         e_1 = y + c ~ breadth,
                                         e_2 = e_1 + c ~ breadth,
                                         e_3 = e_2 + c ~ breadth,
                                         means = e_3 ~ latent_dim,
                                         logvars = e_3 ~ latent_dim,
                                         .act = list(nn_relu,
                                                     logvars = nn_identity,
                                                     means = nn_identity))

                   self$decoder <- nndag(z = ~ latent_dim,
                                         c = ~ c_embed_dim,
                                         d_1 = z + c ~ breadth,
                                         d_2 = d_1 + c ~ breadth,
                                         d_3 = d_2 + c ~ breadth,
                                         out = d_3 ~ input_dim,
                                         .act = list(nn_relu,
                                                     out = nn_identity))
                   
                   self$embedder <- nn_embedding(c_dim, c_embed_dim)

                   self$loggamma <- nn_parameter(torch_tensor(loggamma_init))
                   
                 },
                 reparameterize = function(mean, logvar) {
                   std <- torch_exp(torch_tensor(0.5, device = "cuda") * logvar)
                   eps <- torch_randn_like(std)
                   eps * std + mean
                 },
                 loss_function = function(reconstruction, input, mean, log_var) {
                   kl <- torch_sum(torch_exp(log_var) + torch_square(mean) - log_var, dim = 2L) - latent_dim
                   recon1 <- torch_sum(torch_square(input - reconstruction), dim = 2L) / torch_exp(self$loggamma)
                   recon2 <- self$input_dim * self$loggamma + torch_log(torch_tensor(2 * pi, device = "cuda")) * self$input_dim
                   loss <- torch_mean(recon1 + recon2 + kl)
                   list(loss, torch_mean(recon1*torch_exp(self$loggamma)), torch_mean(kl))
                 },
                 encode = function(x, c = NULL) {
                   if(is.null(c)) {
                     c_embedding <- torch_zeros(x$size()[[1]], self$c_embed_dim, device = x$device)
                   } else {
                     c_embedding <- self$embedder(c)
                   }
                   self$encoder(c_embedding, x)
                 },
                 decode = function(z, c = NULL) {
                   if(is.null(c)) {
                     c_embedding <- torch_zeros(x$size()[[1]], self$c_embed_dim, device = x$device)
                   } else {
                     c_embedding <- self$embedder(c)
                   }
                   self$decoder(c_embedding, z)
                 },
                 forward = function(x, c = NULL) {
                   if(is.null(c)) {
                     c_embedding <- torch_zeros(x$size()[[1]], self$c_embed_dim, device = x$device)
                   } else {
                     c_embedding <- self$embedder(c)
                   }
                   c(means, log_vars) %<-% self$encoder(c_embedding, x)
                   z <- self$reparameterize(means, log_vars)
                   list(self$decoder(c_embedding, z), x, means, log_vars)
                 })
  
  input_dim <- 64L
  c_dim <- 10L
  c_embed_dim <- 64L
  latent_dim <- 64L
  breadth <- 1024L
  
  vae <- vae_mod(input_dim, c_dim, c_embed_dim, latent_dim, breadth, loggamma_init = -3)
  vae <- vae$cuda()
  
  num_epochs <- 25000

  lr <- 0.002
  optimizer <- optim_adam(vae$parameters, lr = lr)
  scheduler <- lr_one_cycle(optimizer, max_lr = lr,
                            epochs = num_epochs, steps_per_epoch = length(train_dl),
                            cycle_momentum = FALSE)

  for (epoch in 1:num_epochs) {

    batchnum <- 0
    coro::loop(for (b in train_dl) {

        batchnum <- batchnum + 1
        optimizer$zero_grad()
        
        c(reconstruction, input, mean, log_var) %<-% vae(b$codes$cuda(), b$troph$cuda())
        c(loss, reconstruction_loss, kl_loss) %<-% vae$loss_function(reconstruction, input, mean, log_var)
        
        ## unconditional model
        #c(reconstruction2, input2, mean2, log_var2) %<-% vae(b$codes$cuda())
        #c(loss2, reconstruction_loss2, kl_loss2) %<-% vae$loss_function(reconstruction2, input2, mean2, log_var2)

        # loss_total <- (loss + loss2) / 2
        # recon_total <- (reconstruction_loss + reconstruction_loss2) / 2
        # kl_total <- (kl_loss + kl_loss2) / 2
        
        if(batchnum %% 2 == 0) {

            cat("Epoch: ", epoch,
                "  batch: ", batchnum,
                "  loss: ", as.numeric(loss$cpu()),
                "  recon loss: ", as.numeric(reconstruction_loss$cpu()),
                "  KL loss: ", as.numeric(kl_loss$cpu()),
                "  loggamma: ", as.numeric(vae$loggamma$cpu()),
                #"    loggamma: ", loggamma,
                "  cond. active dims: ", as.numeric((torch_exp(log_var)$mean(dim = 1L) < 0.5)$sum()$cpu()),
                #"  uncond. active dims: ", as.numeric((torch_exp(log_var2)$mean(dim = 1L) < 0.5)$sum()$cpu()),
                "\n")

        }
        #loss_total$backward()
        loss$backward()
        optimizer$step()
        scheduler$step()
    })
  }
  
  vae

}
