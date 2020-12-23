#' Test for Hardy-Weinberg equilibrium
#'
#' @param x A genind object
#' @param n_boot The number of bootstrap replicates to run
#'
#' @return A tibble with the HWE test results for each population separately
#' @export
#'
#' @importFrom adegenet nPop popNames
#' @importFrom dplyr bind_rows mutate
#' @importFrom pegas hw.test
#' @importFrom tibble rownames_to_column
#'
#' @examples
#' @import adegenet
hwe_test <- function(x, n_boot = 1000) {

  # First, check to see that there are populations defined
  if (is.null(adegenet::pop(x))) {

    # Produce a warning
    warning("Warning: no populations defined. Assuming only one population.")

    # Make a dummy pop "pop1"
    adegenet::pop(x) <- rep("pop1", adegenet::nInd(x))
  }

  # Split the genind object by population

  pop_lst <- adegenet::seppop(x)
  hwe_lst <- purrr::map(pop_lst, function(gen) {

    hw_res <- pegas::hw.test(gen, B = n_boot)
    hw_res <- dplyr::as_tibble(hw_res, rownames = "locus")
    colnames(hw_res) <- c("locus", "chi2", "df", "pval_chi2", "pval_exact")
    hw_res$pop <- as.character(adegenet::pop(gen)[1])
    hw_res <- dplyr::select(hw_res, locus, pop, dplyr::everything())

    return(hw_res)

    })

  # Combine the list of tibbles into a single tibble
  hwe_tbl <- dplyr::bind_rows(hwe_lst)

  return(hwe_tbl)
}
