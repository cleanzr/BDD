#' @include BDDFit.R
NULL

setGeneric("extract", 
           function(x, params = NULL, include = TRUE) standardGeneric("extract"))

#' Extract Samples
#' 
#' @description 
#' Extract samples from a fitted model
#' 
#' @param x an [`BDDFit-class`] object
#' @param params an optional character vector specifying the names of
#'   parameters/summary statistics to extract. If not specified, all 
#'   parameters and summary statistics are extracted.
#' @param include a logical scalar indicating whether the parameters in 
#'   `params` should be included (TRUE) or excluded (FALSE).
#'  
#' @return
#' If params is of length 1, the samples for the requested parameter are 
#' returned as a [`coda::mcmc`] object.
#' 
#' If params is NULL or of length greater than 1, the parameters are returned 
#' in a named list. The samples for each parameter are represented as a
#' [`coda::mcmc`] object.
#' 
#' @examples
#' ## Initialize a BDD model for RLdata500
#' library(comparator) # provides scoring functions
#' 
#' # Add record ids to the dataframe
#' RLdata500$ID <- seq.int(nrow(RLdata500))
#' 
#' # Specify model parameters
#' scoring_fns <- list(
#'   fname_c1 = Levenshtein(normalize=TRUE),
#'   lname_c1 = Levenshtein(normalize=TRUE),
#'   by = BinaryComp(),
#'   bm = BinaryComp(),
#'   bd = BinaryComp()
#' )
#' scoring_breaks <- list(
#'   fname_c1 = c(-Inf, .05, .2, .4, Inf),
#'   lname_c1 = c(-Inf, .05, .2, .4, Inf),
#'   by = c(-Inf, 0, Inf),
#'   bm = c(-Inf, 0, Inf),
#'   bd = c(-Inf, 0, Inf)
#' )
#' lambda <- list(
#'   fname_c1 = c(0.8, 0.85, 0.99),
#'   lname_c1 = c(0.8, 0.85, 0.99),
#'   by = c(0.8, 0.99),
#'   bm = c(0.8, 0.99),
#'   bd = c(0.8, 0.99)
#' )
#' 
#' # Prepare pairwise comparison data
#' pairs <- pairs_all(RLdata500$ID)
#' pairs <- compute_scores(pairs, RLdata500, scoring_fns, id_col='ID')
#' pairs <- discretize_scores(pairs, scoring_breaks)
#' pairs$candidate <- (pairs$fname_c1 < 4) & (pairs$lname_c1 < 4)
#' 
#' model <- BDD(pairs, lambda, candidate_col = 'candidate')
#' 
#' ## Run inference
#' fit <- run_inference(model, n_samples=100, burnin_interval=100)
#' 
#' ## Extract samples of the m probabilities
#' m <- extract(fit, "m")
#' 
#' @aliases extract
#' @export
setMethod("extract", signature = c(x = "BDDFit"), 
          function (x, params = NULL, include = TRUE) {
            if (is.null(params)) return(x@history)
            avail_params <- names(x@history)
            sel_params <- if (include) {
              intersect(params, avail_params)
            } else {
              setdiff(avail_params, params)
            }
            if (length(params) == 1 & length(sel_params) == 1) {
              x@history[[sel_params]]
            } else {
              x@history[sel_params]
            }
          })

#' Extract Complete Samples of the Linkage Structure
#' 
#' @description 
#' This function extracts samples of the linkage structure which cover **all** 
#' records, including records that are excluded from the candidate pairs.
#' 
#' @param x An instance of [`BDDFit-class`]. This is the output 
#'   of [`run_inference`].
#' @param all_rec_ids A vector containing the record identifiers for 
#'   all records.
#' 
#' @return A [`coda::mcmc`] object where the rows correspond to samples and 
#'   the columns correspond to records.
#' 
#' @seealso In contrast to this function, the [`extract`] function extracts 
#'   samples of the linkage structure for a *subset* of records - only those 
#'   records that appear in the candidate pairs.
#' 
#' @examples
#' ## Initialize a BDD model for RLdata500
#' library(comparator) # provides scoring functions
#' 
#' # Add record ids to the dataframe
#' RLdata500$ID <- seq.int(nrow(RLdata500))
#' 
#' # Specify model parameters
#' scoring_fns <- list(
#'   fname_c1 = Levenshtein(normalize=TRUE),
#'   lname_c1 = Levenshtein(normalize=TRUE),
#'   by = BinaryComp(),
#'   bm = BinaryComp(),
#'   bd = BinaryComp()
#' )
#' scoring_breaks <- list(
#'   fname_c1 = c(-Inf, .05, .2, .4, Inf),
#'   lname_c1 = c(-Inf, .05, .2, .4, Inf),
#'   by = c(-Inf, 0, Inf),
#'   bm = c(-Inf, 0, Inf),
#'   bd = c(-Inf, 0, Inf)
#' )
#' lambda <- list(
#'   fname_c1 = c(0.8, 0.85, 0.99),
#'   lname_c1 = c(0.8, 0.85, 0.99),
#'   by = c(0.8, 0.99),
#'   bm = c(0.8, 0.99),
#'   bd = c(0.8, 0.99)
#' )
#' 
#' # Prepare pairwise comparison data
#' pairs <- pairs_all(RLdata500$ID)
#' pairs <- compute_scores(pairs, RLdata500, scoring_fns, id_col='ID')
#' pairs <- discretize_scores(pairs, scoring_breaks)
#' pairs$candidate <- (pairs$fname_c1 < 4) & (pairs$lname_c1 < 4)
#' 
#' model <- BDD(pairs, lambda, candidate_col = 'candidate')
#' 
#' ## Run inference
#' fit <- run_inference(model, n_samples=100, burnin_interval=100)
#' 
#' ## Extract samples of the linkage structure
#' links <- complete_links_samples(fit, RLdata500$ID)
#' 
#' @importFrom coda mcmc
#' @export
complete_links_samples <- function(x, all_rec_ids) {
  if (!inherits(x, 'BDDFit')) 
    stop("`x` must be a BDDFit object")
  if (!is.vector(all_rec_ids))
    stop("`all_rec_ids` must be a vector")
  
  # Original ids of records that were considered as candidate matches
  c_rec_ids <- names(x@state@rec_ids)
  
  # Need all record ids in the same character representation
  all_rec_ids <- as.character(all_rec_ids)
  
  # Check consistency
  missing_rec_ids <- setdiff(c_rec_ids, all_rec_ids)
  if (length(missing_rec_ids) > 0)
    stop("`all_rec_ids` is missing some record ids that appear in `x`")
  
  # Records that were not considered as candidate matches
  nc_rec_ids <- setdiff(all_rec_ids, c_rec_ids)
  
  # Merge
  links <- extract(x, "links")
  num_samples <- nrow(links)
  full_links <- matrix(0L, nrow = num_samples, ncol = length(all_rec_ids), 
                       dimnames = list(NULL, all_rec_ids))
  full_links[,c_rec_ids] <- links
  full_links[,nc_rec_ids] <- matrix(seq_along(nc_rec_ids) + length(c_rec_ids), 
                                    nrow = num_samples, ncol = length(nc_rec_ids), 
                                    byrow = TRUE)
  mcpar <- attr(x@history, 'mcpar')
  full_links <- coda::mcmc(full_links, start=mcpar[1], end=mcpar[2], thin=mcpar[3])
  
  return(full_links)
}