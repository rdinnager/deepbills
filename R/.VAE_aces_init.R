library(fibre)
library(phyf)
library(tidyverse)
library(ape)
library(targets)

bill_z <- read_csv("data/bills_vae_latent_codes_means_only_16dim.csv")
bill_tree <- pf_as_phylo(tar_read(bird_beak_avonet))

bill_pf <- pf_as_pf(bill_tree) %>%
  left_join(bill_z, by = c("label" = "Species"))

ace_fit <- fibre(latent_1 +
                   latent_2 +
                   latent_3 + 
                   latent_4 +
                   latent_5 +
                   latent_6 +
                   latent_7 +
                   latent_8 +
                   latent_9 +
                   latent_10 +
                   latent_11 +
                   latent_12 +
                   latent_13 +
                   latent_14 +
                   latent_15 +
                   latent_16 ~
                   bre_brownian(phlo),
                 data = bill_pf,
                 verbose = 2,
                 engine_options = list(control.family = list(hyper = list(hyper = list(prec = list(prior = "pc.prec", initial = 4, fixed = TRUE))))))

rates <- ace_fit$random$phlo %>%
  select(ID, mean) %>%
  separate(ID, c("var", "node"), ":") %>%
  pivot_wider(names_from = var, values_from = mean) %>%
  select(node,
         latent_1,
         latent_2,
         latent_3,
         latent_4,
         latent_5,
         latent_6,
         latent_7,
         latent_8,
         latent_9,
         latent_10,
         latent_11,
         latent_12,
         latent_13,
         latent_14,
         latent_15,
         latent_16)

ace_preds <- predict(ace_fit)
names(ace_preds) <- paste0(names(ace_preds), "_pred")
ace_means <-  map(ace_preds, ".pred_mean") %>%
  as_tibble()

bill_w_preds <- bill_pf %>%
  bind_cols(ace_means) %>%
  mutate(time = pf_flow_sum(phlo))

write_rds(bill_w_preds, "data/bill_w_preds_16dim.rds")

bill_edge_trajs <- pf_ends(bill_w_preds$phlo) %>%
  mutate(isna = is.na(end),
         end = ifelse(isna, start, end),
         start = ifelse(isna, NA, start)) %>%
  left_join(bill_w_preds %>%
              select(start = label,
                     ends_with("_pred")) %>%
              rename_with(~ paste0(.x, "_start"),
                          ends_with("_pred"))) %>%
  left_join(bill_w_preds %>%
              select(end = label,
                     ends_with("_pred")) %>%
              rename_with(~ paste0(.x, "_end"),
                          ends_with("_pred"))) %>%
  left_join(bill_w_preds %>%
              select(start = label,
                     start_time = time)) %>%
  left_join(bill_w_preds %>%
              select(end = label,
                     end_time = time)) %>%
  mutate(start_time = ifelse(isna, 0, start_time),
         across(ends_with("_start"), ~ ifelse(is.na(.x), 0, .x)))

write_rds(bill_edge_trajs, "data/bill_edge_trajs_16dim.rds")

ggplot(bill_edge_trajs,
       aes(x = latent_1_pred_start, y = latent_2_pred_start)) +
  geom_segment(aes(xend = latent_1_pred_end, yend = latent_2_pred_end)) +
  coord_equal() +
  theme_minimal() +
  geom_point(aes(latent_1_pred, latent_2_pred), data = bill_w_preds %>% select(-phlo) %>% as_tibble(),
             colour = "red", size = 0.1)

ggplot(bill_edge_trajs,
       aes(x = latent_15_pred_start, y = latent_16_pred_start)) +
  geom_segment(aes(xend = latent_15_pred_end, yend = latent_16_pred_end)) +
  coord_equal() +
  theme_minimal() +
  geom_point(aes(latent_15_pred, latent_16_pred), data = bill_w_preds %>% dplyr::filter(is_tip) %>% select(-phlo) %>% as_tibble(),
             colour = "red", size = 0.25)

bill_rates <- bill_edge_trajs %>%
  mutate(vecs = pick(ends_with("_end")) - pick(ends_with("_start")),
         time_len = end_time - start_time) %>%
  mutate(vecs = vecs / time_len) %>%
  select(label = end, vecs) %>%
  unnest(vecs)

colnames(bill_rates) <- gsub("_pred_end", "", colnames(bill_rates))

write_rds(bill_pf, "data/bill_pf_16dim.rds")
write_csv(bill_rates, "data/init_rates_16dim.csv")