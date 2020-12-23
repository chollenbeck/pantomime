#' Run a quick and dirty PCA using a genind object
#'
#' @param x A genind oject
#' @param n_axes The number of axes to keep
#'
#' @return A tidy data frame with individual loadings
#' @export
#'
#' @importFrom adegenet nInd pop scaleGen
#' @importFrom ade4 dudi.pca
#' @importFrom dplyr as_tibble everything mutate select
#'
#' @examples
qpca <- function(x, n_axes = 4) {

  if (is.null(adegenet::pop(x))) {
    pop(gen) <- rep("pop1", adegenet::nInd(x))
  }

  # Scale the data using the mean
  x_gen <- adegenet::scaleGen(x, NA.method = "mean")

  # Run the PCA with adegenet
  pca_obj <- ade4::dudi.pca(x_gen, center = FALSE, scale = FALSE, scannf = FALSE, nf = n_axes)

  # Convert the PCA object to a data frame
  pca_tbl <- as.data.frame(pca_obj$li)
  pca_tbl <- cbind(pca_tbl, pop = as.character(adegenet::pop(x)))

  # Add the individual names from the PCA
  pca_tbl <- dplyr::as_tibble(pca_tbl, rownames = "ind")
  pca_tbl <- dplyr::mutate(pca_tbl, pop = as.character(pop))

  pca_tbl <- dplyr::select(pca_tbl, .data$ind, .data$pop, dplyr::everything())

  return(pca_tbl)

}
