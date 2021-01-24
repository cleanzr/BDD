#' @include utils.R
NULL

.check_BDDModel <- function(object) {
  # TODO: implement
  # level_counts, agree_levels must be a list of integer vectors, should have 
  # same dimensions and same names
  # pair_index, row_id_index should be a list of vectors
  # nc_counts should contain "alpha" and "beta" sublists, each of which has 
  # integer vectors for each attribute.
  # priors should contain "alpha0", "alpha1", "beta0", "beta1", "lambda" 
  # sublists, each with strictly positive numeric vectors for each attribute.
  # candidate_pairs should contain two columns: "ID_1" and "ID_2". Should be 
  # factors?
  # comparisons should contain a column for each attribute. Must be integer.
  # m, u: list with entry for each attribute. Must be a numeric vector. 
  # length is one less than number of score levels.
  # links: length must be same as length of rec_ids. Values between 0 
  # and ..?
  # iteration: non-negative scalar.
  return(TRUE)
}

#' Container class for the model
#' 
#' @description
#' This class records the state of model parameters that are realized 
#' during MCMC sampling. It also stores indexing data structures for 
#' used to accelerate sampling. 
#' 
#' @slot iteration Non-negative integer. Counts the number of successive 
#'   applications of the Markov transition operator. This is not necessarily 
#'   equal to the number of samples if thinning is applied.
#' @slot m Named list containing the `m*` probabilities for each attribute. 
#'   Each entry corresponds to an attribute in `comparisons` and contains a 
#'   vector of `m*` probabilities for the discrete score levels. Due to the 
#'   sequential parameterization of `m*`, there is no probability included 
#'   for the largest score level.
#' @slot u a named list containing the `u*` probabilities for each attribute. 
#'   Each entry corresponds to an attribute in `comparisons` and contains a 
#'   vector of `u*` probabilities for the discrete score levels. Due to the 
#'   sequential parameterization of `u*`, there is no probability included 
#'   for the largest score level.
#' @slot comparisons a data frame of discrete comparison scores. Each row 
#'   corresponds to a pair of records. The discrete comparison scores for 
#'   each attribute are stored in columns whose names match the names in 
#'   `m`, `u`, `level_counts`, etc. 
#' @slot candidate_pairs a data frame that specifies pairs of records that 
#'   are candidate matches. It contains two columns, `ID_1` and `ID_2`, which 
#'   contain the record identifiers for each pair. Each pair is included 
#'   only once: e.g. if pair "A, B" appears, then "B, A" does not appear.
#' @slot priors a list containing the priors for `alpha1`, `alpha0`, 
#'   `beta1`, `beta0` and `lambda`.
#' @slot nc_counts a list containing summary statistics for record pairs that 
#'   are non-candidate matches. These are used when updating the `u*` 
#'   probabilities.
#' @slot row_id_index An index from rec_ids to row indices in `candidate_pairs` 
#'   where they appear. This is used to speed up inference.
#' @slot pair_index An index for `candidate_pairs`. Maps a record identifier 
#'   to a vector of record identifiers that are potential matches (candidate
#'   pairs).
#' @slot level_counts A named list containing counts for each attribute 
#'   at each level of agreement for candidate matches.
#' @slot rec_ids A named vector that maps the original record identifiers 
#'   (the names) to standardized integer identifiers (the entries). Records 
#'   that are not candidate matches are excluded.
#' @slot links An integer vector of entity identifiers for each record in 
#'   `rec_ids`. The i-th entry gives the linked entity for the 
#'   i-th record in `rec_ids`.
#' 
#' @seealso 
#' The [`BDD`] function should be used to initialize a model 
#' from data.
BDDModel <- setClass("BDDModel", slots = c(iteration = "integer",
                                           links = "integer",
                                           m = "list",
                                           u = "list",
                                           comparisons = "data.frame",
                                           candidate_pairs = "data.frame",
                                           priors = "list",
                                           nc_counts = "list",
                                           row_id_index = "list",
                                           pair_index = "list",
                                           level_counts = "list",
                                           rec_ids = "vector"),
                     validity = .check_BDDModel)

#' @importFrom utils capture.output
setMethod("show", "BDDModel", function(object) {
  n_c_pairs <- nrow(object@candidate_pairs) # number of candidate pairs
  n_nc_pairs <- (object@nc_counts$A0[[1]][1] + 
                   object@nc_counts$cum_A0[[1]][1]) # number of non-candidate pairs
  n_records <- length(object@rec_ids)
  level_counts <- object@level_counts
  cat("BDDModel\n",
      "Defined on ", n_c_pairs + n_nc_pairs, " record pairs: ", n_c_pairs, 
      " candidate and ", n_nc_pairs, " non-candidate matches.\n",
      "Attributes used for matching:\n", sep="")
  for (a in names(level_counts)) {
    cat("  * ", a, ": with ", length(level_counts[[a]]), " levels of agreement\n", sep="")
  }
})

#' Initialize a Bayesian Duplicate Detection Model
#' 
#' Initializes a [`BDDModel-class`] given observed data and model 
#' hyperparameters.
#' 
#' @param comparisons A data frame containing discrete comparison scores for 
#'   record pairs. Each row corresponds to a record pair whose identifiers 
#'   are specified in the columns named in `id_cols`. The discrete comparison 
#'   scores for each attribute are stored in columns whose names match the 
#'   names in `lambda`, `alpha1`, `beta1`, etc. An optional column (named in 
#'   `candidate_col`) may be included, which specifies whether the record 
#'   pair is a candidate match or not. If this column is omitted, all record 
#'   pairs are considered to be candidate matches. Any additional columns are 
#'   ignored. 
#' @param lambda A named list containing the lower truncation points for the 
#'   truncated Beta priors on the `m*` probabilities. Each element of the list 
#'   corresponds to an attribute in `comparisons` and must share the same
#'   name. The value of each list element must be a numeric vector of 
#'   probabilities (on the open unit interval), whose i-th entry corresponds 
#'   to the i-th level of agreement for that attribute. Due to the sequential 
#'   parameterization of m*, a prior is not required for the highest level of 
#'   agreement.
#' @param alpha1 A named list containing "alpha" shape parameters for the 
#'   truncated Beta priors on the `m*` probabilities. The specifications of 
#'   the list are the same as for `lambda`, except the numeric vectors must be 
#'   strictly positive. If NULL, defaults to 1 for each attribute/level of 
#'   agreement.
#' @param beta1 A named list containing "beta" shape parameters for the 
#'   truncated Beta priors on the `m*` probabilities. The specifications of 
#'   the list are the same as for `lambda`, except the numeric vectors must be 
#'   strictly positive. If NULL, defaults to 1 for each attribute/level of 
#'   agreement.
#' @param alpha0 A named list containing "alpha" shape parameters for the 
#'   Beta priors on the `u*` probabilities. The specification is the same as 
#'   for `alpha1`. If NULL, defaults to 1 for each attribute/level of 
#'   agreement.
#' @param beta0 A named list containing "beta" shape parameters for the 
#'   Beta priors on the `u*` probabilities. The specification is the same as 
#'   for `beta1`. If NULL, defaults to 1 for each attribute/level of 
#'   agreement.
#' @param id_cols A character vector specifying the pair of column names in 
#'   `comparisons` that contain the record identifiers. Defaults to 
#'   `c("ID.x", "ID.y")`.
#' @param candidate_col An optional string referencing a column in 
#'   `comparisons` that encodes whether each pair is a candidate match
#'   or not. Defaults to NULL.
#' 
#' @return a [`BDDModel-class`] object. This can be used as the starting 
#' point for inference using the [`run_inference`] function.
#' 
#' @references
#' Sadinle, Mauricio. 2014. "Detecting Duplicates in a Homicide Registry 
#' Using a Bayesian Partitioning Approach." *Ann. Appl. Stat.* **8** (4): 
#' 2404–34. <https://doi.org/10.1214/14-AOAS779>.
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
#' model <- BDD(pairs, lambda)
#' 
#' @importFrom stats na.omit
#' @importFrom utils head relist
#' @export
BDD <- function(comparisons, lambda, alpha1 = NULL, beta1 = NULL, 
                alpha0 = NULL,  beta0 = NULL, id_cols = c("ID.x", "ID.y"), 
                candidate_col = NULL) {
  if (!is.character(id_cols) || length(id_cols) != 2) 
    stop("`id_cols` must be a character vector of length 2")
  
  attributes <- names(lambda)
  req_colnames <- c(id_cols, attributes, candidate_col)
  missing_colnames <- setdiff(req_colnames, colnames(comparisons))
  if (length(missing_colnames) > 0)
    stop("`comparisons` is missing the following columns: ", paste0(missing_colnames, collapse = ", "))
  
  init_shape_param <- function(shape, name) {
    if (!is.null(shape)) {
      if (!all(unlist(shape) > 0)) 
        stop(paste0("`", name, "` must contain numeric vectors with strictly positive values"))
    } else {
      # Set to default
      shape <- lapply(lambda, function(x) rep(1, times = length(x)))
    }
    shape
  }
  alpha0 <- init_shape_param(alpha0, "alpha0")
  alpha1 <- init_shape_param(alpha1, "alpha1")
  beta0 <- init_shape_param(beta0, "beta0")
  beta1 <- init_shape_param(beta1, "beta1")
  
  lambda_flat <- unlist(lambda)
  if (!all(lambda_flat >= 0 & lambda_flat <= 1))
    stop("`lambda` must contain numeric vectors with values on the unit interval [0, 1]")
  
  # Consistency checks and extracting agreement levels
  agree_levels <- lapply(lambda, function(x) seq_len(length(x) + 1))
  
  for (attribute in attributes) {
    n_agree_lvls <- length(lambda[[attribute]]) + 1
    
    if ( (n_agree_lvls - 1 != length(alpha1[[attribute]])) || 
         (n_agree_lvls - 1!= length(beta1[[attribute]])) ) {
      stop("mismatch in dimensions of alpha1, beta1 and lambda")
    }
    
    domain <- na.omit(unique(comparisons[[attribute]]))
    invalid_vals <- setdiff(domain, seq_len(n_agree_lvls))
    if (length(invalid_vals) > 0) {
      stop("attribute '", attribute, "' contains invalid values")
    }
  }
  
  # Deal with prior non-matches (non-candidate pairs)
  if (!is.null(candidate_col)) {
    if (!is.character(candidate_col) || length(candidate_col) != 1) {
      stop("`candidate_col` must be a character vector of length 1")
    }
    
    candidate <- comparisons[[candidate_col]]
    # Counts at each level of agreement for *non-candidate* pairs (which are 
    # guaranteed non-matches)
    A0 <- lapply(attributes, function(a) {
      tabulate(comparisons[!candidate, a], nbins = length(lambda[[a]]) + 1)
    })
    names(A0) <- attributes
    
    cum_A0 <- lapply(A0, function(x) rev(cumsum(rev(x))))
    A0 <- lapply(A0, function(x) x[-length(x)])
    cum_A0 <- lapply(cum_A0, function(x) x[-1])
    
    # Remove non-candidate pairs to save memory (not needed beyond this point)
    comparisons <- comparisons[candidate,]
  } else {
    A0 <- lapply(lambda, function(x) rep(0L, times=length(x)))
    cum_A0 <- lapply(lambda, function(x) rep(0L, times=length(x)))
  }
  
  # Copy rec_ids and attributes into separate dataframes
  candidate_pairs <- subset(comparisons, select=id_cols)
  colnames(candidate_pairs) <- c("ID_1", "ID_2")
  comparisons <- comparisons[attributes]
  
  # Get set of original rec_ids that appear in comparisons
  orig_rec_ids <- unique(unlist(candidate_pairs))
  orig_rec_ids <- sort(orig_rec_ids)
  
  # Use new (more convenient) rec_ids that start at 1
  rec_ids <- seq_along(orig_rec_ids)
  names(rec_ids) <- as.character(orig_rec_ids)

  # Map original rec_ids to new rec_ids
  candidate_pairs <- lapply(candidate_pairs, function(x) {
    unclass(factor(x, orig_rec_ids))
  })
  candidate_pairs <- as.data.frame(candidate_pairs)
  
  # Build an index for the candidate pairs. Maps a record id to a vector of 
  # rec_ids that are candidate pairs.
  pair_index <- lapply(rec_ids, function(rec_id) {
    pair_ids <- c(candidate_pairs$ID_1[candidate_pairs$ID_2 == rec_id], 
                  candidate_pairs$ID_2[candidate_pairs$ID_1 == rec_id])
    pair_ids <- pair_ids[pair_ids != rec_id]
    # Transform to index starting at zero for C++ code
    return(sort(pair_ids) - 1L)
  })
  
  # Build an index from rec_ids to row ids in `candidate_pairs` where they 
  # appear
  row_id_index <- lapply(rec_ids, function(rec_id) {
    # Transform to index starting at zero for C++ code
    which( (candidate_pairs$ID_1==rec_id)|(candidate_pairs$ID_2==rec_id) ) - 1L
  })

  # For each attribute, get counts at each level of agreement for all 
  # record pairs
  level_counts <- lapply(attributes, function(x) {
    tabulate(comparisons[[x]], nbins = length(agree_levels[[x]]))
  })
  names(level_counts) <- attributes

  # Initialize unobserved variables
  m <- rep(NA_real_, times=sum(sapply(lambda, length)))
  m <- relist(m, lambda)
  links <- seq.int(from = 0L, length.out = length(rec_ids))
  names(links) <- names(rec_ids)
  
  model <- BDDModel(iteration = 0L,
                    links = links,
                    m = m,
                    u = m,
                    comparisons = comparisons,
                    candidate_pairs = candidate_pairs,
                    priors = list(alpha0=alpha0, alpha1=alpha1, 
                                  beta0=beta0, beta1=beta1, lambda=lambda),
                    nc_counts = list(A0=A0, cum_A0=cum_A0),
                    row_id_index = row_id_index,
                    pair_index = pair_index,
                    level_counts = level_counts,
                    rec_ids = rec_ids)
  
  return(model)
}