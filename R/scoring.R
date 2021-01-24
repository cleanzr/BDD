
#' Compute Pairwise Comparison Scores
#' 
#' Given a set of record pairs, this function computes comparison scores 
#' for each attribute. This is known as the "scoring" stage of a record 
#' linkage/entity resolution pipeline.
#' 
#' @param pairs A data frame or matrix specifying pairs of records in `x` to 
#'   compare, where rows index pairs and columns index the identifiers of the 
#'   two constituent records. If `id_col` is specified, the record are 
#'   identified according to the corresponding column in `x`. Otherwise if 
#'   `id_col=NULL` the records are identifier by their row index.
#' @param x A data frame or matrix of records. Columns must be present for 
#'   each attribute referenced in `scoring_fns`. An optional column containing 
#'   record identifiers is used if `id_col` is not NULL. Any other columns 
#'   are ignored.
#' @param scoring_fns A list containing scoring functions for the attributes 
#'   to compare. Each element of the list must be named according to an 
#'   attribute in `x`. The scoring function must take a pair of vectors as its 
#'   first two arguments, and return a vector of scores.
#' @param id_col An integer or string specifying a column in `x` that contains 
#'   record identifiers. Defaults to NULL.
#' 
#' @return The `pairs` data frame augmented with columns containing scores 
#'   for each of the attributes in `scoring_fns`.
#' 
#' @seealso 
#' The [`discretize_scores`] function can be used to discretize the scores 
#' data frame output by this function.
#' To prepare the record `pairs` data frame, refer to the indexing 
#' functions: [`pairs_all`], [`pairs_hamming`] or [`pairs_fuzzyblock`].
#' 
#' @importFrom progress progress_bar
#' @export
compute_scores <- function(pairs, x, scoring_fns, id_col = NULL) {
  # Input validation
  if (ncol(pairs) != 2) stop("`pairs` must contain exactly two columns")
  
  if (!is.list(scoring_fns) && !all(sapply(scoring_fns, is.function))) {
    stop("`scoring_fns` must be a list of functions")
  }
  if (!is.null(id_col) & length(id_col) != 1) {
    stop("`id_col` must be a vector of length 1 if not NULL")
  }
  
  attributes <- names(scoring_fns)
  req_colnames <- attributes
  if (!is.null(id_col) & is.character(id_col)) {
    req_colnames <- c(req_colnames, id_col)
  }
  missing_colnames <- setdiff(attributes, colnames(x))
  if (length(missing_colnames) > 0) {
    stop("`x` is missing the following columns: ", paste0(missing_colnames, collapse = ", "))
  }
  
  # Initialize data frame for scores
  scored_pairs <- as.data.frame(pairs)
  
  # Use id_col to name columns if not NULL
  if (!is.null(id_col)) {
    colnames(scored_pairs) <- paste0(id_col, c(".x", ".y"))
  }
  
  # Display progress
  pb <- progress::progress_bar$new(total = length(attributes), clear = TRUE)
  pb$tick(0)
  for (attribute in attributes) {
    id_left <- scored_pairs[[1]]
    id_right <- scored_pairs[[2]]
    if (!is.null(id_col)) {
      # Convert to row index if using id_col
      id_left <- match(id_left, x[[id_col]])
      id_right <- match(id_right, x[[id_col]])
    }
    values_left <- x[[attribute]][id_left]
    values_right <- x[[attribute]][id_right]
    scored_pairs[[attribute]] <- scoring_fns[[attribute]](values_left, values_right)
    pb$tick()
  }
  pb$terminate()
  
  scored_pairs
}

#' Discretize Pairwise Comparison Scores
#' 
#' Given a set of scored record pairs, this function converts the scores 
#' to discrete levels.
#' 
#' @param pairs A data frame of scored record pairs as returned by 
#'   [`compute_scores`]. Must contain the attributes referenced in 
#'   `breaks`.
#' @param breaks A list of scoring breaks for each scores attribute in 
#'   `pairs`. The names of the entries must match the colnames in `pairs`. 
#'   The scoring breaks should be represented as a vector of breaks or cut 
#'   points, which specify intervals closed on the right.
#' 
#' @return A transformed version of the input `pairs` data frame, in which 
#'   scored columns are discretized. Note that the discrete scores are 
#'   encoded as positive integers, where 1 represents the smallest score. 
#' 
#' @export
discretize_scores <- function(pairs, breaks) {
  if (!is.list(breaks) && !all(sapply(breaks, is.numeric))) {
    stop("`breaks` must be a list of numeric vectors")
  }
  
  attributes <- names(breaks)
  missing_colnames <- setdiff(attributes, colnames(pairs))
  if (length(missing_colnames) > 0) {
    stop("`pairs` is missing the following columns: ", paste0(missing_colnames, collapse = ", "))
  }
  
  for (attribute in attributes) {
    pairs[[attribute]] <- cut(pairs[[attribute]], breaks = breaks[[attribute]], labels = FALSE)
  }
  
  return(pairs)
}