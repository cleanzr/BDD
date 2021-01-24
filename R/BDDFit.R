#' @include BDDModel.R utils.R
NULL

#' Fitted Model
#' 
#' @description 
#' Stores data associated with a fitted [`BDDModel-class`] object, 
#' including posterior samples, diagnostics, etc.
#' 
#' @slot history A list containing model parameters of interest and summary 
#'   statistics along the Markov chain. See below for details.
#' @slot state A [`BDDModel-class`] object which represents the state 
#'   of the model at the end of the Markov chain.
#' 
#' @details 
#' The `history` may include the following [`coda::mcmc`] objects:
#' \describe{
#'   \item{links}{A matrix recording samples of the linkage/clustering 
#'     structure, where rows correspond to samples and columns correspond to 
#'     records. For efficiency reasons, links are only included for 
#'     records in the set of candidate pairs. Records not in the set of 
#'     candidate pairs are constrained to occupy their own cluster (i.e. they 
#'     are prevented from linking with other records). Each sample of the 
#'     linkage/clustering structure is encoded as an integer membership 
#'     vector: records with the same integer value are assigned to the same 
#'     cluster.}
#'   \item{m}{A matrix object recording the m*-probabilities for each 
#'     attribute/agreement level. Rows index samples along the Markov chain 
#'     and columns index attributes/agreement levels.}
#'   \item{u}{A matrix object recording the u*-probabilities for each 
#'     attribute/agreement level. Rows index samples along the Markov chain 
#'     and columns index attributes/agreement levels.}
#'   \item{n_clusters}{A vector object recording the number of 
#'     clusters/entities the records have been partitioned into. This number 
#'     only accounts for records in the set of candidate pairs---i.e. it 
#'     excludes singleton clusters associated with records in the set of 
#'     non-candidate pairs.}
#' }
#' They can be accessed using the [`extract`] method.
setClass("BDDFit", 
         slots = c(history = "list",
                   state = "BDDModel"))

#' @importFrom mcmcse multiESS ess
setMethod("show", "BDDFit", function(object) {
  # Extract info about burn-in, thinning, number of samples
  mcpar <- attr(object@history, "mcpar")
  start_iter <- mcpar[1]
  end_iter <- mcpar[2]
  thin <- mcpar[3]
  n_samples <- end_iter - start_iter + thin
  
  # Compute effective sample size for particular variables, if present
  valid_ess_varnames <- c("n_linked_ents", "m", "u")
  ess <- list()
  for (varname in intersect(valid_ess_varnames, names(object@history))) {
    var <- object@history[[varname]]
    if (ncol(var) > 1)
      ess[[varname]] <- tryCatch(mcmcse::multiESS(var), error = function(e) NaN)
    else
      ess[[varname]] <- tryCatch(mcmcse::ess(var), error = function(e) NaN)
  }
  
  cat("Fitted BDDModel\n",
      n_samples, if (n_samples > 1) " samples" else " sample", 
      " after thinning (interval=", thin, 
      ") with burn-in (interval=", start_iter, ")\n",
      "Recorded parameters/summary statistics: \n", 
      "  ",  toString(names(object@history)), "\n", sep="")
  if (length(ess) > 0) {
    cat("Effective sample size:\n")
    for (varname in names(ess)) {
      cat("  ", varname, ": ", ess[[varname]], "\n", sep="")
    }
  }
})

#' Function to concatenate two BDDFit objects
#' 
#' @param resultA,resultB [`BDDFit-class`] objects. `resultB` must occur 
#'   after `resultA` along the Markov chain. Additionally, the results must use 
#'   the same `thin_interval`.
#' @return an [`BDDFit-class`] object
#' 
#' @importFrom coda mcmc
#' @noRd
combine_results <- function(resultA, resultB) {
  # Assume that both inputs are instances of type 'BDDFit' and that 
  # they refer to the same data set & model
  
  mcparA <- attr(resultA@history, 'mcpar')
  mcparB <- attr(resultB@history, 'mcpar')
  
  # Check validity of inputs
  if (diff(mcparA[1:2]) == 0) stop("resultA has no history")
  if (diff(mcparB[1:2]) == 0) stop("resultB has no history")
  if (mcparA[2] >= mcparB[1] ) {
    stop("resultB must occur chronologically after resultA")
  }
  if (mcparA[3] != mcparB[3]) 
    stop("results must have the same thin_interval")
  
  if ( !setequal(names(resultA@history), names(resultB@history)) ) {
    stop("results must have the same history variables")
  }
  
  # Combine history matrices/vectors
  hist_combined <- vector(mode="list", length = length(resultA@history))
  attr(hist_combined, 'mcpar') <- c(mcparA[1], resultB@state@iteration, mcparA[3])
  names(hist_combined) <- names(resultA@history)
  for (varname in names(resultA@history)) {
    varA <- resultA@history[[varname]]
    varB <- resultB@history[[varname]]
    if (is.matrix(varA)) {
      history <- rbind(varA, varB)
    } else {
      history <- c(varA, varB)
    }
    hist_combined[[varname]] <- coda::mcmc(history, start=mcparA[1], thin = mcparA[3])
  }
  
  return(new("BDDFit", history=hist_combined, state=resultB@state))
}