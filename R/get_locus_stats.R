#' Calculate summary statistics by locus for each population in a genind object
#'
#' @param x A genind object
#' @param hwe_test Calculate p-value for HWE
#' @param n_boot Number of bootstrap replicates for HWE test
#'
#' @return A tibble with summary statistics (by locus and population)
#' @export
#'
#' @importFrom adegenet genind2genpop makefreq nInd pop seppop
#' @importFrom dplyr as_tibble bind_rows case_when group_by left_join mutate rename select slice ungroup
#' @importFrom purrr map map2
#' @importFrom tibble data_frame rownames_to_column
#' @importFrom tidyr pivot_longer separate
#'
#' @examples
get_locus_stats <- function(x, hwe_test = FALSE, n_boot = 1000) {

  #browser()

  # First, check to see that there are populations defined
  if (is.null(adegenet::pop(x))) {

    # Produce a warning
    warning("Warning: no populations defined. Assuming only one population.")

    # Make a dummy pop "pop1"
    adegenet::pop(x) <- rep("pop1", adegenet::nInd(x))
  }

  # Make a tibble of minor allele frequency
  maf_tbl <- adegenet::genind2genpop(x, quiet = TRUE)
  maf_tbl <- adegenet::makefreq(maf_tbl, quiet = TRUE)
  maf_tbl <- tibble::as_tibble(maf_tbl, rownames = "pop")
  maf_tbl <- tidyr::pivot_longer(maf_tbl, cols = -c("pop"), names_to = "locus", values_to = "freq")
  maf_tbl <- tidyr::separate(maf_tbl, locus, into = c("locus", "allele"), sep = "\\.")
  maf_tbl <- dplyr::group_by(maf_tbl, pop, locus)
  maf_tbl <- dplyr::slice(maf_tbl, which.min(freq))
  maf_tbl <- dplyr::select(maf_tbl, -allele)
  maf_tbl <- dplyr::rename(maf_tbl, maf = freq)
  maf_tbl <- dplyr::ungroup(maf_tbl)

  # Make a tibble for allele counts
  allele_tbl <- adegenet::genind2genpop(x, quiet = TRUE)
  allele_tbl <- adegenet::makefreq(allele_tbl, quiet = TRUE)
  allele_tbl <- dplyr::as_tibble(allele_tbl, rownames = "pop")
  allele_tbl <- tidyr::pivot_longer(allele_tbl, -("pop"), names_to = "locus", values_to = "freq")
  allele_tbl <- tidyr::separate(allele_tbl, locus, into = c("locus", "allele"), sep = "\\.")
  allele_tbl <- dplyr::filter(allele_tbl, freq > 0)
  allele_tbl <- dplyr::group_by(allele_tbl, pop)
  allele_tbl <- dplyr::count(allele_tbl, locus)
  allele_tbl <- dplyr::select(allele_tbl, locus, pop, n_alleles = n)
  allele_tbl <- dplyr::ungroup(allele_tbl)

  # Create a tibble for missing data
  pop_lst <- adegenet::seppop(x)
  loc_typed <- purrr::map(pop_lst, ~adegenet::propTyped(., by = "loc"))
  prop_typed <- purrr::map(loc_typed, ~tibble::tibble(locus = names(.), prop_typed = .))
  md_lst <- purrr::map2(prop_typed, names(prop_typed), ~dplyr::mutate(.x, pop = .y))
  md_tbl <- dplyr::bind_rows(md_lst)
  md_tbl <- dplyr::mutate(md_tbl, prop_missing = 1 - prop_typed)
  md_tbl <- dplyr::select(md_tbl, locus, pop, prop_missing)

  # Create a tibble for heterozygosity
  het_lst <- purrr::map(pop_lst, function(gen) {

      hf <- hierfstat::genind2hierfstat(gen)

      # Fix the problem where hierfstat renames the locus when there is only one locus
      if (ncol(hf) == 2) {
        colnames(hf)[2] <- locNames(gen)
      }

      stats <- hierfstat::basic.stats(hf)
      stat_tbl <- stats$perloc
      stat_tbl <- dplyr::as_tibble(stat_tbl, rownames = "locus")

      stat_tbl$pop <- as.character(adegenet::pop(gen)[1])
      stat_tbl <- dplyr::select(stat_tbl, locus, pop, ho = Ho, he = Hs, fis = Fis)

      return(stat_tbl)

  })

  # Combine the list of tibbles into a single tibble
  het_tbl <- dplyr::bind_rows(het_lst)

  # Combine the data into a tibble
  fin_tbl <- dplyr::left_join(maf_tbl, allele_tbl, by = c("locus", "pop"))
  fin_tbl <- dplyr::left_join(fin_tbl, md_tbl, by = c("locus", "pop"))
  fin_tbl <- dplyr::left_join(fin_tbl, het_tbl, by = c("locus", "pop"))
  fin_tbl <- dplyr::select(fin_tbl, locus, pop, n_alleles, prop_missing, maf, ho, he, fis)
  fin_tbl <- dplyr::mutate(fin_tbl, maf = dplyr::case_when(n_alleles == 1 ~ 0,
                                                           n_alleles == 2 ~ maf,
                                                           n_alleles > 2 ~ maf))

  if (hwe_test == TRUE) {
    hwe_tbl <- hwe_test(x, n_boot = n_boot)
    hwe_tbl <- dplyr::select(hwe_tbl, locus, pop, hwe_pval = pval_exact)

    fin_tbl <- dplyr::left_join(fin_tbl, hwe_tbl, by = c("locus", "pop"))
  } else {
    fin_tbl$hwe_pval <- NA
  }

  return(fin_tbl)

}
