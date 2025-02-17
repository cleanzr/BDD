% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BDDModel.R
\name{BDD}
\alias{BDD}
\title{Initialize a Bayesian Duplicate Detection Model}
\usage{
BDD(
  comparisons,
  lambda,
  alpha1 = NULL,
  beta1 = NULL,
  alpha0 = NULL,
  beta0 = NULL,
  id_cols = c("ID.x", "ID.y"),
  candidate_col = NULL
)
}
\arguments{
\item{comparisons}{A data frame containing discrete comparison scores for
record pairs. Each row corresponds to a record pair whose identifiers
are specified in the columns named in \code{id_cols}. The discrete comparison
scores for each attribute are stored in columns whose names match the
names in \code{lambda}, \code{alpha1}, \code{beta1}, etc. An optional column (named in
\code{candidate_col}) may be included, which specifies whether the record
pair is a candidate match or not. If this column is omitted, all record
pairs are considered to be candidate matches. Any additional columns are
ignored.}

\item{lambda}{A named list containing the lower truncation points for the
truncated Beta priors on the \verb{m*} probabilities. Each element of the list
corresponds to an attribute in \code{comparisons} and must share the same
name. The value of each list element must be a numeric vector of
probabilities (on the open unit interval), whose i-th entry corresponds
to the i-th level of agreement for that attribute. Due to the sequential
parameterization of m*, a prior is not required for the highest level of
agreement.}

\item{alpha1}{A named list containing "alpha" shape parameters for the
truncated Beta priors on the \verb{m*} probabilities. The specifications of
the list are the same as for \code{lambda}, except the numeric vectors must be
strictly positive. If NULL, defaults to 1 for each attribute/level of
agreement.}

\item{beta1}{A named list containing "beta" shape parameters for the
truncated Beta priors on the \verb{m*} probabilities. The specifications of
the list are the same as for \code{lambda}, except the numeric vectors must be
strictly positive. If NULL, defaults to 1 for each attribute/level of
agreement.}

\item{alpha0}{A named list containing "alpha" shape parameters for the
Beta priors on the \verb{u*} probabilities. The specification is the same as
for \code{alpha1}. If NULL, defaults to 1 for each attribute/level of
agreement.}

\item{beta0}{A named list containing "beta" shape parameters for the
Beta priors on the \verb{u*} probabilities. The specification is the same as
for \code{beta1}. If NULL, defaults to 1 for each attribute/level of
agreement.}

\item{id_cols}{A character vector specifying the pair of column names in
\code{comparisons} that contain the record identifiers. Defaults to
\code{c("ID.x", "ID.y")}.}

\item{candidate_col}{An optional string referencing a column in
\code{comparisons} that encodes whether each pair is a candidate match
or not. Defaults to NULL.}
}
\value{
a \code{\linkS4class{BDDModel}} object. This can be used as the starting
point for inference using the \code{\link{run_inference}} function.
}
\description{
Initializes a \code{\linkS4class{BDDModel}} given observed data and model
hyperparameters.
}
\examples{
## Initialize a BDD model for RLdata500
library(comparator) # provides scoring functions

# Add record ids to the dataframe
RLdata500$ID <- seq.int(nrow(RLdata500))

# Specify model parameters
scoring_fns <- list(
  fname_c1 = Levenshtein(normalize=TRUE),
  lname_c1 = Levenshtein(normalize=TRUE),
  by = BinaryComp(),
  bm = BinaryComp(),
  bd = BinaryComp()
)
scoring_breaks <- list(
  fname_c1 = c(-Inf, .05, .2, .4, Inf),
  lname_c1 = c(-Inf, .05, .2, .4, Inf),
  by = c(-Inf, 0, Inf),
  bm = c(-Inf, 0, Inf),
  bd = c(-Inf, 0, Inf)
)
lambda <- list(
  fname_c1 = c(0.8, 0.85, 0.99),
  lname_c1 = c(0.8, 0.85, 0.99),
  by = c(0.8, 0.99),
  bm = c(0.8, 0.99),
  bd = c(0.8, 0.99)
)

# Prepare pairwise comparison data
pairs <- pairs_all(RLdata500$ID)
pairs <- compute_scores(pairs, RLdata500, scoring_fns, id_col='ID')
pairs <- discretize_scores(pairs, scoring_breaks)
pairs$candidate <- (pairs$fname_c1 < 4) & (pairs$lname_c1 < 4)

model <- BDD(pairs, lambda)

}
\references{
Sadinle, Mauricio. 2014. "Detecting Duplicates in a Homicide Registry
Using a Bayesian Partitioning Approach." \emph{Ann. Appl. Stat.} \strong{8} (4):
2404–34. \url{https://doi.org/10.1214/14-AOAS779}.
}
