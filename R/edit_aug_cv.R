#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title
#' @param trophic_niche_cv_ai64_aug
#' @return
#' @author rdinnager
#' @export
edit_aug_cv <- function(trophic_niche_cv_ai64_aug) {

  #split <- trophic_niche_cv_ai64_aug[[1]][[1]]
  add_out_ids <- function(split) {
    diffs <- setdiff(1:nrow(split$data), split$in_id)
    split$out_id <- diffs[which(split$data$type[diffs] == "original")]
    split
  }
  trophic_niche_cv_ai64_aug[[1]] <- map(trophic_niche_cv_ai64_aug[[1]], add_out_ids)

  trophic_niche_cv_ai64_aug
  
}
