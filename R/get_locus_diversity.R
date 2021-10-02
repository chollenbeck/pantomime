#' Generate a tidy data frame of diversity stats for each locus
#'
#' @param x A genind object
#'
#' @return A tibble
#'
#' Columns are:
#' n_alleles: the number of observed alleles in the population
#' shannon: shannon diversity (calculated from poppr)
#' hexp_nei: Nei's gene diversity
#' evenness: a measure of allelic evenness, with one being completely even representation (calculated from poppr)
#'
#'
#' @export
#'
#' @importFrom adegenet nInd pop seppop
#' @importFrom dplyr bind_rows if_else filter mutate select
#' @importFrom poppr locus_table
#' @importFrom purrr map
#' @importFrom rlang .data
#' @importFrom tibble as_tibble
#'
#' @examples
get_locus_diversity <- function(x) {


  # First, check to see that there are populations defined
  if (is.null(adegenet::pop(x))) {

    # Produce a warning
    warning("Warning: no populations defined. Assuming only one population.")

    # Make a dummy pop "pop1"
    adegenet::pop(x) <- rep("pop1", adegenet::nInd(x))
  }

  pops <- adegenet::seppop(x)

  div_lst <- purrr::map(pops, function(i) {
      loc_tbl <- tibble::as_tibble(poppr::locus_table(i, index = "shannon", information = FALSE), rownames = "locus")

      loc_tbl <- dplyr::mutate(loc_tbl, pop = popNames(i))
      loc_tbl <- dplyr::filter(loc_tbl, .data$locus != "mean")
      loc_tbl <- dplyr::mutate(loc_tbl, n_alleles = as.integer(.data$allele),
                                        shannon = as.numeric(.data$H),
                                        hexp_nei = as.numeric(.data$Hexp),
                                        evenness = as.numeric(.data$Evenness))

      loc_tbl <- dplyr::select(loc_tbl, .data$locus, .data$pop, .data$n_alleles, .data$shannon, .data$hexp_nei, .data$evenness)

      loc_tbl <- dplyr::mutate(loc_tbl, evenness = dplyr::if_else(is.infinite(.data$evenness), NA_real_, .data$evenness))


    })

  return(dplyr::bind_rows(div_lst))

}
