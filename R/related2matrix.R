#' Convert output from related 'coancestry' fuction into a matrix
#'
#' @param x Output from related's coancestry() function
#' @param method Relatedness estimator to use
#' @param rows A character vector to subset the rows of the matrix
#' @param cols A character vector to subset the columns of the matrix
#'
#' @return A matrix of relatedness estimates
#' @export
#'
#' @importFrom dplyr mutate select
#' @importFrom rlang .data enquo
#' @importFrom tidyr pivot_wider
#' @importFrom tidyselect all_of
#'
#'
#' @examples
related2matrix <- function(x,
                              method = c("ritland", "wang", "lynchli", "lynchrd", "quellergt", "dyadml", "trioml"),
                              rows = NULL,
                              cols = NULL) {


  # Check the parameters

  # Check the row input
  if (!is.null(rows)) {
    if (!is.character(rows)) {
      stop("Row individuals need to be provided in the form of a character vector")
    }
  }

  # Check the column input
  if (!is.null(cols)) {
    if (!is.character(cols)) {
      stop("Column individuals need to be provided in the form of a character vector")
    }
  }

  # Select the appropriate columns of the related output
  met <- rlang::enquo(method)
  relate_tbl <- dplyr::select(x$relatedness, .data$ind1.id, .data$ind2.id, !!met)

  # Define the correct order of individuals
  indiv_order <- unique(c(relate_tbl$ind1.id, relate_tbl$ind2.id))

  # Convert the table into a pairwise matrix
  rel_mat <- tidyr::pivot_wider(relate_tbl, names_from = .data$ind2.id, values_from = !!met)
  rel_mat <- dplyr::mutate(rel_mat, !!indiv_order[1] := NA)
  rel_mat <- dplyr::select(rel_mat, .data$ind1.id, tidyselect::all_of(indiv_order))

  rel_mat <- rel_mat[match(indiv_order, rel_mat$ind1.id),]
  rel_mat <- dplyr::select(rel_mat, -.data$ind1.id)
  rownames(rel_mat) <- indiv_order

  rel_mat[lower.tri(rel_mat)] <- t(rel_mat)[lower.tri(rel_mat)]

  # Check if row and column subsets are defined and apply them
  if (!is.null(rows)) {
    rwz <- rows
  } else {
    rwz <- indiv_order
  }

  if (!is.null(cols)) {
    clz <- cols
  } else {
    clz <- indiv_order
  }

  return(rel_mat[rwz, clz])

}
