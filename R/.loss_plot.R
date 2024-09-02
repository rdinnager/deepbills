library(tidyverse)

losses <- read_delim("data/model/sdf_net_training.csv",
                     col_names = c("Epoch", "Duration", "Mean Loss", "Latent Variance")) |>
  mutate(epoch = 25000 - (n():1))

plot(losses$`Mean Loss`)

ggplot(losses, aes(x = epoch, y = `Mean Loss`)) +
  geom_line() +
  scale_x_continuous(labels = scales::label_comma()) +
  labs(x = "Epoch",
       y = "Mean Loss") +
  theme_minimal()