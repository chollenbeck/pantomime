#' Convert a genind object to 'related' format
#'
#' @param x A genind object
#'
#' @return List object used as input to 'related'
#' @export
#'
#' @importFrom adegenet genind2df indNames
#' @importFrom dplyr across bind_cols funs select transmute_all
#' @importFrom related readgenotypedata
#' @importFrom rlang .data
#' @importFrom tidyselect everything
#'
#' @examples
  genind2related <- function(x) {

  # Convert the genind object to a dataframe
  rel <- dplyr::bind_cols(inds = adegenet::indNames(x), adegenet::genind2df(x, oneColPerAll = TRUE))

  # Related needs missing data coded as '0'. If one of the alleles is already 0, it needs to be recoded:
  # Check the genind object to see if there are alleles coded as zero
  if ("0" %in% unlist(alleles(x))) {
    rel <-  dplyr::mutate(rel, dplyr::across(.cols = c(-inds, -pop), .fns = ~as.integer(.x) + 1))
    rel <- dplyr::mutate(rel, dplyr::across(.cols = tidyselect::everything(), .fns = ~gsub('NA', 0, .x)))
    rel <- related::readgenotypedata(dplyr::select(rel, -pop))
  } else {
    rel <- dplyr::mutate(rel, across(.cols = tidyselect::everything(), .fns = ~gsub('NA', 0, .x)))
    rel <- related::readgenotypedata(dplyr::select(rel, -pop))
  }

  return(rel)
}
