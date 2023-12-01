#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_conf_mat2
make_conf_mat_plot <- function(trophic_niche_conf_mat2) {

  conf_df <- trophic_niche_conf_mat2$table %>%
    as.data.frame() %>%
    mutate(correct = factor(ifelse(Prediction == Truth, "correct", "incorrect"),
                            levels = c("incorrect", "correct")))
  
  ggplot(conf_df,
       aes(y = Freq, axis1 = Truth, axis2 = Prediction)) +
    geom_alluvium(aes(fill = correct), width = 1/2) +
    geom_stratum(width = 1/12, fill = "grey80", color = "black") +
    #geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
    geom_label_repel(stat = "stratum", aes(label = after_stat(stratum)),
                     min.segment.length = 0, point.padding = NA, 
                     label.padding = 0.25,
                     fill = scales::alpha("white", 0.4)) +
    scale_x_discrete(limits = c("Truth", "Prediction"), expand = c(0, 0)) +
    scale_fill_brewer(type = "qual", palette = "Set1", name = "") +
    theme_minimal() +
    theme(axis.text = element_text(size = 12))

}
