#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param div_rate_clads
#' @param bird_beak_codes
reconcile_clads <- function(div_rate_clads, bird_beak_codes) {

  tre <- pf_as_phylo(div_rate_clads)
  tre <- drop.tip(tre, tre$tip.label[!tre$tip.label %in% bird_beak_codes$label],
                  trim.internal = TRUE, collapse.singles = FALSE,
                  subtree = FALSE)
  ord <- rev(postorder(tre))
  #tre <- reorder(tre, "postorder")
  
  edge_degree <- ape::degree(tre, details = TRUE)[tre$edge[ , 2]][ord]
  drle <- rle(edge_degree)
  srle <- drle
  srle$values <- as.character(srle$values)
  srle$values[drle$values == 2] <- paste0("sr_", seq_along(drle$values[srle$values == 2]))
  and_back <- inverse.rle(srle)
  and_back[grep("sr_", and_back) + 1] <- and_back[grep("sr_", and_back)]
  and_back[and_back == "1"] <- paste0("tr_", seq_along(and_back[and_back == "1"]))
  and_back[and_back == "3"] <- paste0("br_", seq_along(and_back[and_back == "3"]))
  
  singleton_df <- tibble(cat = and_back, 
                         node_name = c(tre$tip.label, tre$node.label)[tre$edge[ , 2]][ord],
                         len = tre$edge.length[ord]) %>%
    left_join(div_rate_clads %>%
                select(-phlo, -is_tip, node_name = label, div_rate)) %>%
    group_by(cat) %>%
    mutate(new_div_rate = weighted.mean(div_rate, len))

  tre1 <- pf_as_phylo(div_rate_clads)
  tre1 <- drop.tip(tre1, tre1$tip.label[!tre1$tip.label %in% bird_beak_codes$label])
  
  tre2 <- tre1
  tre2 <- makeNodeLabel(tre2)
  matcher <- tibble(orig = c(tre1$tip.label, tre1$node.label), 
                    new = c(tre2$tip.label, tre2$node.label)) %>%
    left_join(singleton_df %>%
                ungroup() %>%
                select(orig = node_name, new_div_rate))
  # %>%
  #   left_join(div_rate_clads %>%
  #               select(-phlo, -is_tip, orig = label, div_rate))
  tre_pf <- pf_as_pf(tre1) %>%
    left_join(matcher, by = c("label" = "new")) %>%
    left_join(bird_beak_codes %>%
                select(label, starts_with("latent_code_"))) %>%
    left_join(div_rate_clads %>%
                select(label, tip_rate))
  
  tre_pf

}
