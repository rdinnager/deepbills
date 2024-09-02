library(tidyverse)
library(targets)

tar_load(ai_codes)
tar_load(bad_birds)
tar_load(trophic_niche_dat_ai_all)

ai_codes <- ai_codes |>
  filter(!Species %in% bad_birds)

hist(ai_codes$latent_1, breaks = 100, main = "Latent 1", xlab = "Value")

pdf("output/ai_codes_gaussianicity.pdf", width = 10, height = 10)

oldpar <- par(mfrow = c(4, 4))

for (i in 1:64) {
  hist(trophic_niche_dat_ai64_all[[paste0("latent_code_", i)]], breaks = 100, main = paste("DeepSDF Latent", i), xlab = "Value")
}

for (i in 1:64) {
  qqnorm(trophic_niche_dat_ai64_all[[paste0("latent_code_", i)]], main = paste("DeepSDF Latent", i))
  qqline(trophic_niche_dat_ai64_all[[paste0("latent_code_", i)]])
}

par(mfrow = c(4, 4))

for (i in 1:15) {
  hist(ai_codes[[paste0("latent_", i)]], breaks = 100, main = paste("VAE Stage2 Latent", i), xlab = "Value")
}

plot.new()

for (i in 1:15) {
  qqnorm(ai_codes[[paste0("latent_", i)]], main = paste("VAE Stage2 Latent", i))
  qqline(ai_codes[[paste0("latent_", i)]])
}

dev.off()