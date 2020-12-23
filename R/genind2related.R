#' Convert a genind object to 'related' format
#'
#' @param x A genind object
#'
#' @return List object used as input to 'related'
#' @export
#'
#' @importFrom adegenet genind2df indNames
#' @importFrom dplyr bind_cols funs select transmute_all
#' @importFrom related readgenotypedata
#' @importFrom rlang .data
#'
#' @examples
genind2related <- function(x) {

  rel <- dplyr::bind_cols(inds = adegenet::indNames(x), adegenet::genind2df(x, oneColPerAll = TRUE))
  rel <- dplyr::transmute_all(rel, dplyr::funs(gsub('NA', 0, .data$.)))
  rel <- related::readgenotypedata(dplyr::select(rel, -pop))
  return(rel)
}
