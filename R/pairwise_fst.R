#' Calculate Weir and Cockerham's Fst between pairs of populations
#'
#' @param x A genind object
#' @param n_bootstraps Number of bootstraps to calculate confidence intervals
#'
#' @return A dataframe with pairwise FST, lower and upper 95% confidence limits
#' @export
#'
#' @importFrom dplyr filter left_join rename
#' @importFrom hierfstat boot.ppfst genind2hierfstat pairwise.WCfst
#' @importFrom rlang .data
#' @importFrom tibble as_tibble
#' @importFrom tidyr pivot_longer
#'
#' @examples
pairwise_fst <- function(x, n_bootstraps = 1000) {

  hf <- hierfstat::genind2hierfstat(x)

  fst <- hierfstat::pairwise.WCfst(hf, diploid = TRUE)

  fst_tbl <- tibble::as_tibble(fst, rownames = "pop")
  fst_tbl <- tidyr::pivot_longer(fst_tbl, -.data$pop, names_to = "pop2", values_to = "fst")
  fst_tbl <- dplyr::rename(fst_tbl, pop1 = .data$pop)

  fst_cf <- hierfstat::boot.ppfst(hf, nboot = n_bootstraps)

  fst_ll <- tibble::as_tibble(fst_cf$ll, rownames = "pop")
  fst_ll <- tidyr::pivot_longer(fst_ll, -.data$pop, names_to = "pop2", values_to = "ll")
  fst_ll <- dplyr::rename(fst_ll, pop1 = .data$pop)
  fst_ll <- dplyr::filter(fst_ll, !is.na(.data$ll))

  fst_ul <- tibble::as_tibble(fst_cf$ul, rownames = "pop")
  fst_ul <- tidyr::pivot_longer(fst_ul, -.data$pop, names_to = "pop2", values_to = "ul")
  fst_ul <- dplyr::rename(fst_ul, pop1 = .data$pop)
  fst_ul <- dplyr::filter(fst_ul, !is.na(.data$ul))

  fst_tbl <- dplyr::left_join(fst_tbl, fst_ll, by = c("pop1", "pop2"))
  fst_tbl <- dplyr::left_join(fst_tbl, fst_ul, by = c("pop1", "pop2"))

  return(fst_tbl)

}
