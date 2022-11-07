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

  stat_tbl <- pop_lst %>%
    map(function(i) {
        fis <- hierfstat::genind2hierfstat(i) %>%
          as.data.frame() %>%
          tibble::rownames_to_column(var = "id") %>%
          #tibble::as.tibble() %>%
          dplyr::filter(id != "dumind") %>%
          dplyr::select(-pop) %>%
          hierfstat::betas(betaijT = TRUE) %>%
          .$betaij %>%
          diag() * 2 - 1

        # Make a data frame of genotypes
        gt_tbl <- x@tab %>%
          tibble::as_tibble(rownames = "ind") %>%
          tidyr::pivot_longer(-ind, names_to = "locus_allele", values_to = "gt") %>%
          tidyr::extract(locus_allele, c("locus", "allele"), "(.*)\\.(.*)")

        # Make a heterozygosity table
        het_tbl <- gt_tbl %>%
          dplyr::mutate(het = dplyr::if_else(gt == 1, 1, 0)) %>%
          dplyr::group_by(ind) %>%
          dplyr::summarize(hz = sum(het, na.rm = TRUE) / dplyr::n())

        comb_tbl <- dplyr::left_join(tibble::tibble(ind = names(fis), fis = fis), het_tbl, by = "ind")

        return(comb_tbl)

    }) %>%
    dplyr::bind_rows()

  ind_tbl <- dplyr::left_join(typed_tbl, stat_tbl, by = "ind")


  return(ind_tbl)

}
