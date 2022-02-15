#' Calculate summary statistics by individual using a genind object
#'
#' @param x A genind object
#'
#' @return A tibble of summary stats for individuals
#' @export
#'
#' @importFrom adegenet nInd pop propTyped seppop
#' @importFrom dplyr bind_rows filter left_join mutate select
#' @importFrom hierfstat betas genind2hierfstat
#' @importFrom purrr map map2
#' @importFrom tibble rownames_to_column tibble
#'
#' @examples
get_ind_stats <- function(x) {

  # First, check to see that there are populations defined
  if (is.null(adegenet::pop(x))) {

    # Produce a warning
    warning("Warning: no populations defined. Assuming only one population.")

    # Make a dummy pop "pop1"
    adegenet::pop(x) <- rep("pop1", adegenet::nInd(x))
  }

  # Get a missing data tibble
  pop_lst <- adegenet::seppop(x)
  typed_lst <- purrr::map(pop_lst, ~adegenet::propTyped(., by = "ind"))
  typed_lst <- purrr::map(typed_lst, ~tibble::tibble(ind = names(.), prop_typed = .))
  typed_lst <- purrr::map2(typed_lst, names(typed_lst), ~dplyr::mutate(.x, pop = .y))
  typed_tbl <- dplyr::bind_rows(typed_lst)
  typed_tbl <- dplyr::mutate(typed_tbl, prop_missing = 1 - prop_typed)
  typed_tbl <- dplyr::select(typed_tbl, ind, pop, prop_missing)

  fis_tbl <- pop_lst %>%
    map(function(i) {
        fis <- hierfstat::genind2hierfstat(i) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "id") %>%
          dplyr::filter(id != "dumind") %>%
          dplyr::select(-pop) %>%
          hierfstat::betas(betaijT = TRUE) %>%
          .$betaij %>%
          diag() * 2 - 1

          return(tibble::tibble(ind = names(fis), fis = fis))

    }) %>%
    dplyr::bind_rows()

  ind_tbl <- dplyr::left_join(typed_tbl, fis_tbl, by = "ind")


  return(ind_tbl)

}
