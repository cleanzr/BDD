#' @include utils.R BDDFit.R
NULL

#' Run Inference for a Model
#' 
#' @description
#' Runs approximate inference for a model using Markov chain Monte 
#' Carlo (MCMC).
#' 
#' @param x an \R object representing a model or a model with sampling 
#'   history. Currently, there are methods defined for [`BDDModel-class`] 
#'   objects and [`BDDFit-class`] objects. Passing a model begins 
#'   sampling from scratch, in which case a burn-in period is recommended. 
#'   Passing a model with sampling history will resume inference, adding 
#'   samples to the existing history.
#' @param n_samples a positive integer. The number of samples to generate 
#'   after applying burn-in and thinning.
#' @param thin_interval a positive integer. The period for saving samples 
#'   (intermediate samples are discarded). The default value is 1, which 
#'   means no thinning is applied.
#' @param burnin_interval a non-negative integer. The number of initial samples 
#'   to discard as burn-in. The default value is 0, which means no burn-in is 
#'   applied.
#' @param ... further arguments passed to or from other methods. 
#' @return 
#' Returns a [`BDDFit-class`] object with the following slots:
#' \item{history}{a list containing the sampling history for 
#'   parameters/summary statistics along the Markov chain.}
#' \item{state}{a model object of the same class as the input, which 
#'   represents the state of all parameters in the last step of the Markov 
#'   chain.}
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
#' @rdname run_inference
#' @export
setGeneric("run_inference", 
           function(x, n_samples, thin_interval=1, burnin_interval=0, ...) {
             standardGeneric("run_inference")
           })

#' @describeIn run_inference Specialization for [`BDDModel-class`]
#' @importFrom coda mcmc
#' @export
setMethod("run_inference", signature(x = "BDDModel"), 
  function(x, n_samples, thin_interval, burnin_interval, ...) {
    # Check validity of arguments
    if (!is.scalar(n_samples)) stop("n_samples must be a scalar")
    if (!is.scalar(thin_interval)) stop("thin_interval must be a scalar")
    if (!is.scalar(burnin_interval)) stop("burnin_interval must be a scalar")
    if (n_samples < 0) stop("n_samples must be positive")
    if (thin_interval < 1) stop("thin_interval must be positive")
    if (burnin_interval < 0) stop("burnin_interval must be non-negative")
    if (as.integer(n_samples) != n_samples) 
      stop("n_samples must be an integer")
    if (as.integer(thin_interval) != thin_interval) 
      stop("thin_interval must be an integer")
    if (as.integer(burnin_interval) != burnin_interval) 
      stop("burnin_interval must be an integer")
    
    t1 <- Sys.time()
    if (x@iteration == 0) {
      # Not resuming a Markov chain
      start_iter <- ifelse(burnin_interval, burnin_interval, thin_interval)
    } else {
      # Resuming a Markov chain
      start_iter <- x@iteration + thin_interval
    }
    
    # Call C++ function
    result <- .sample(x, n_samples, thin_interval, burnin_interval)
    
    # These attributes are needed for the resuming logic
    attr(result@history, 'mcpar') <- c(start_iter, result@state@iteration, thin_interval)
    
    # Convert history to coda::mcmc objects
    for (i in seq_along(result@history)) {
      result@history[[i]] <- coda::mcmc(result@history[[i]], start=start_iter, 
                                        thin=thin_interval)
    }
    t2 <- Sys.time()
    message("Completed sampling in ", format(t2-t1))
    
    return(result)
  })

#' @describeIn run_inference Specialization for [`BDDFit-class`]
#' @export
setMethod("run_inference", signature(x = "BDDFit"), 
  function(x, n_samples, thin_interval, burnin_interval, ...) {
    # Resuming a previous chain
    if (burnin_interval != 0) 
      warning("ignoring burnin_interval=", burnin_interval, " when resuming")
    old_thin_interval <- attr(x@history, 'mcpar')[3]
    if (thin_interval != old_thin_interval) 
      warning("using thin_interval=", old_thin_interval, " from previous run")
    result <- run_inference(x@state, n_samples, thin_interval, burnin_interval)
    result <- combine_results(x, result)
    return(result)
  })