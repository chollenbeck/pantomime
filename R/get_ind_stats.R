#' Calculate summary statistics by individual using a genind object
#'
#' @param x A genind object
#'
#' @return A tibble of summary stats for individuals
#' @export
#'
#' @importFrom adegenet nInd pop propTyped seppop
#' @importFrom dplyr bind_rows mutate select
#' @importFrom purrr map map2
#' @importFrom tibble tibble
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




  return(typed_tbl)

}
