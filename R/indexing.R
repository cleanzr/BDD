#' Index All Record Pairs
#' 
#' Given a collection of \eqn{n} records, this function builds an index 
#' containing all \eqn{n(n-1)/2} record pairs. 
#' 
#' @param rec_ids A vector of record identifiers.
#' 
#' @return A data frame of record pairs, where rows index pairs and columns 
#'   index the identifiers of records in each pair.
#' 
#' @seealso 
#' [`pairs_hamming`] and [`pairs_fuzzyblock`] are alternative indexing 
#' functions.
#' 
#' @export
pairs_all <- function(rec_ids) {
  x <- rep(list(rec_ids), 2)
  names(x) <- c("ID.x", "ID.y")
  
  # Generate all n^2 pairs
  pairs <- expand.grid(x, stringsAsFactors = FALSE, KEEP.OUT.ATTRS = FALSE)
  
  # Discard duplicate pairs to leave n(n-1)/2 pairs
  pairs[pairs[[1]] < pairs[[2]], ]
}


#' Index Record Pairs Based on Hamming Distance
#' 
#' Given a collection of records, this function builds an index containing 
#' all pairs whose attribute-level Hamming distance falls below some cut-off. 
#' For instance, with a cut-off of 2, all pairs are included that disagree 
#' on 2 or fewer attributes.
#' 
#' @param x A data frame of records. Columns must be present for 
#'   each attribute referenced in `attribute_cols`. An optional column 
#'   containing record identifiers is used if `id_col` is not NULL. Any other 
#'   columns are ignored.
#' @param attr_cols A character vector of attribute columns in `x` to be 
#'   compared when computing the Hamming distance.
#' @param dist_cutoff An integer in the interval `[0, length(attribute_cols)]`. 
#'   All record pairs with Hamming distance less than or equal to this cutoff 
#'   are included in the index.
#' @param id_col A string specifying a column in `x` that contains record 
#'   identifiers. Defaults to NULL.
#' 
#' @return A data frame of record pairs, where rows index pairs and columns 
#'   index the identifiers of records in each pair.
#' 
#' @seealso 
#' [`pairs_all`] and [`pairs_fuzzyblock`] are alternative 
#' indexing functions
#' 
#' @importFrom fuzzyjoin fuzzy_join
#' @export
pairs_hamming <- function(x, attr_cols, dist_cutoff, id_col = NULL) {
  # Validate input
  if (!is.character(attr_cols) || length(attr_cols) < 1) 
    stop("`attr_cols` must be a character vector of length > 0")
  if (!is.null(id_col) & length(id_col) != 1) 
    stop("`id_col` must be a vector of length 1 if not NULL")
  if (!as.integer(dist_cutoff) == dist_cutoff || 
      dist_cutoff > length(attr_cols) || dist_cutoff < 0) {
    stop("`dist_cutoff` must be an integer on the interval [0, ", length(attr_cols), "]")
  }
  
  req_colnames <- c(attr_cols, id_col)
  missing_colnames <- setdiff(req_colnames, colnames(x))
  if (length(missing_colnames) > 0) {
    stop("`x` is missing the following columns: ", paste0(missing_colnames, collapse = ", "))
  }
  
  # Add integer record ids if none provided
  if (is.null(id_col)) {
    id_col <- "ID"
    x[[id_col]] <- seq.int(nrow(x))
  }
  
  multi_match_fun <- function(x, y) rowSums(x != y, na.rm = TRUE) <= dist_cutoff
  multi_by <- attr_cols
  names(multi_by) <- attr_cols
  pairs <- fuzzyjoin::fuzzy_join(x, x, multi_by = multi_by, 
                                 multi_match_fun = multi_match_fun,
                                 mode = "inner")
  
  # Keep only record ids
  pairs <- pairs[paste0(id_col, c(".x", ".y"))]
  
  # Discard duplicate pairs
  pairs[pairs[[1]] < pairs[[2]], ]
}


#' Index Record Pairs Using Fuzzy Blocking
#' 
#' Prepares an index of record pairs that are "similar" as determined by a 
#' distance function applied to one of the attributes. This can be viewed as 
#' a fuzzy generalization of standard blocking.
#' 
#' @param x A data frame of records. Columns must be present for 
#'   each attribute referenced in `attribute_cols`. An optional column 
#'   containing record identifiers is used if `id_col` is not NULL. Any other 
#'   columns are ignored.
#' @param block_col A string specifying a column in `x` to use as the blocking 
#'   variable.
#' @param dist_fn A vectorized distance function.
#' @param dist_cutoff A positive distance cutoff. All record pairs with 
#'   a distance less than or equal to this cutoff are included in the index.
#' @param id_col A string specifying a column in `x` that contains record 
#'   identifiers. Defaults to NULL.
#' @param include_na A logical specifying whether to include pairs for which 
#'   the distance function evaluates to NA. Defaults to FALSE.
#' 
#' @return A data frame of record pairs, where rows index pairs and columns 
#'   index the identifiers of records in each pair.
#' 
#' @seealso 
#' [`pairs_all`] and [`pairs_hamming`] are alternative indexing 
#' functions
#' 
#' @importFrom fuzzyjoin fuzzy_join
#' @export
pairs_fuzzyblock <- function(x, block_col, dist_fn, dist_cutoff, id_col = NULL,
                             include_na = FALSE) {
  # Check input validity
  if (length(block_col) != 1) stop("`block_col` must be a vector of length 1")
  if (!is.null(id_col) & length(id_col) != 1 | !is.character(id_col)) 
    stop("`id_col` must be a character vector of length 1 if not NULL")
  if (!is.function(dist_fn)) stop("`dist_fn` must be a function")
  if (!is.numeric.scalar(dist_cutoff) || dist_cutoff <= 0) 
    stop("`dist_cutoff` must be a positive numeric")
  
  req_colnames <- block_col
  if (!is.null(id_col)) req_colnames <- c(req_colnames, id_col)
  missing_colnames <- setdiff(req_colnames, colnames(x))
  if (length(missing_colnames) > 0) {
    stop("`x` is missing the following columns: ", paste0(missing_colnames, collapse = ", "))
  }
  
  # Add integer record ids if none provided
  if (is.null(id_col)) {
    id_col <- "ID"
    x[[id_col]] <- seq.int(nrow(x))
  }
  
  # Discard unneeded columns
  x <- x[c(id_col, block_col)]
  
  match_fun <- function(x, y) dist_fn(x,y) <= dist_cutoff
  by <- block_col
  names(by) <- block_col
  pairs <- fuzzyjoin::fuzzy_join(x, x, by = by, 
                                 match_fun = match_fun,
                                 mode = "inner")
  
  # Keep only record ids
  pairs <- pairs[paste0(id_col, c(".x", ".y"))]
  
  # Discard duplicate pairs
  pairs[pairs[[1]] < pairs[[2]], ]
}