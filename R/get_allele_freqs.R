#' Get allele frequencies from a genind object
#'
#' @param x A genind object
#'
#' @return A tidy tibble of allele frequencies
#' @export
#'
#' @importFrom adegenet genind2genpop tab
#' @importFrom tidyr extract pivot_longer
#' @importFrom tibble rownames_to_column
#'
#' @examples
get_allele_freqs <- function(x) {

  freq_tbl <- adegenet::genind2genpop(x, quiet = TRUE)
  freq_tbl <- adegenet::tab(freq_tbl, freq = TRUE)
  freq_tbl <- as.data.frame(freq_tbl)
  freq_tbl <- tibble::rownames_to_column(freq_tbl, var = "pop")
  freq_tbl <- tidyr::pivot_longer(freq_tbl, names_to = "allele", values_to = "freq", -.data$pop)
  freq_tbl <- tidyr::extract(freq_tbl, .data$allele, c("locus", "allele"), "(.*)\\.(.*)")

  return(freq_tbl)

}
