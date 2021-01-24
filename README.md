<!-- README.md is generated from README.Rmd. Please edit that file -->

# BDD: Duplicate detection using a Bayesian partitioning approach

This R package implements a Bayesian model for duplicate detection as
proposed in (Sadinle 2014). It resembles the classic Fellegi-Sunter
model (Fellegi and Sunter 1969) in that coreference decisions are based
on pairwise attribute-level comparisons. However, unlike the
Fellegi-Sunter mdoel, this model targets a *partition* of the records,
rather than pairwise coreference predictions, thereby ensuring
transitivity. Inference is conducted using a Gibbs sampler, with parts
implemented in [Rcpp](http://www.rcpp.org/).

## Installation

`BDD` is not currently available on CRAN. The latest development version
can be installed from source using `devtools` as follows:

``` r
library(devtools)
install_github("cleanzr/BDD")
```

## Example

We demonstrate how to perform duplicate detection using the `RLdata500`
test data set included with the package. `RLdata500` contains 500
synthetic personal records, 50 of which are duplicates with
randomly-generated errors. In the code block, below we load the package
and examine the first few rows of `RLdata500`.

``` r
library(magrittr)    # pipe operator (must be loaded before BDD)
library(comparator)  # normalized Levenshtein distance
library(clevr)       # evaluation functions
library(BDD)

RLdata500[['rec_id']] <- seq.int(nrow(RLdata500))
head(RLdata500)
#>   fname_c1 fname_c2 lname_c1 lname_c2   by bm bd rec_id
#> 1  CARSTEN     <NA>    MEIER     <NA> 1949  7 22      1
#> 2     GERD     <NA>    BAUER     <NA> 1968  7 27      2
#> 3   ROBERT     <NA> HARTMANN     <NA> 1930  4 30      3
#> 4   STEFAN     <NA>    WOLFF     <NA> 1957  9  2      4
#> 5     RALF     <NA>  KRUEGER     <NA> 1966  1 13      5
#> 6  JUERGEN     <NA>   FRANKE     <NA> 1929  7  4      6
```

The model uses agreement levels between attributes in order to assess
whether a pair of records is a match (referring to the same entity) or
not. To compute the agreement levels, we’ll first compute real-valued
scores between the attributes, then maps the scores to discrete levels
of agreement. In the code block below, we specify scoring functions for
the five attributes we will use for matching (we don’t use `fname_c2`
and `lname_c2` as more than 90% of the values are missing). Note that
the scoring functions must accept a pair of attribute vectors `x` and
`y` as arguments and return a vector of scores.

``` r
scoring_fns <- list(
  fname_c1 = Levenshtein(normalize = TRUE),
  lname_c1 = Levenshtein(normalize = TRUE),
  by = function(x, y) abs(x - y),
  bm = function(x, y) abs(x - y),
  bd = function(x, y) abs(x - y)
)
```

For each scoring function above, we provide a breaks vector which
specifies the discrete levels of agreement (from ‘high’ agreement to
‘low’).

``` r
scoring_breaks <- list(
  fname_c1 = c(-Inf, .05, .2, .4, Inf),
  lname_c1 = c(-Inf, .05, .2, .4, Inf),
  by = c(-Inf, 0, 1, 3, Inf),
  bm = c(-Inf, 0, 1, 3, Inf),
  bd = c(-Inf, 0, 2, 7, Inf)
)
```

Now we are ready to compute the agreement levels for the record pairs.
Since this is a small data set, we consider all pairs using the
`pairs_all` function. For larger data sets, blocking/indexing is
recommended using `pairs_hamming`, `pairs_fuzzyblock` or a custom
indexing function.

``` r
pairs <- pairs_all(RLdata500$rec_id) %>%
  compute_scores(RLdata500, scoring_fns, id_col = 'rec_id') %>% 
  discretize_scores(scoring_breaks)
```

To speed up inference, we only consider a subset of the pairs as
candidate matches. Specifically, we consider pairs that have a strong
agreement on name (accounting for missing names).

``` r
pairs[['candidate']] <- (pairs$fname_c1 < 4) & (pairs$lname_c1 < 4) | 
                            is.na(pairs$fname_c1) | is.na(pairs$lname_c1) 
```

Next we specify the priors on the m\* and u\* probabilities for each
attribute and agreement level. `lambda` specifies the lower truncation
points for the truncated Beta priors on the m\* probabilities. By
default, a uniform prior is used over the truncated interval; a
non-uniform prior can be specified using the `alpha1` and `beta1`
arguments to `BDD` below. The Beta priors on the u\* probabilities are
also uniform by default and can be adjusted using the `alpha0` and
`beta0` arguments.

``` r
lambda <- list(
  fname_c1 = c(0.8, 0.85, 0.99),
  lname_c1 = c(0.8, 0.85, 0.99),
  by = c(0.8, 0.85, 0.99),
  bm = c(0.8, 0.85, 0.99),
  bd = c(0.8, 0.85, 0.99)
)
```

Finally, we initialize the model and run inference using Markov chain
Monte Carlo.

``` r
model <- BDD(pairs, lambda, id_cols = c("rec_id.x", "rec_id.y"), 
             candidate_col = "candidate")
fit <- run_inference(model, 100, thin_interval = 10, burnin_interval = 100)
#> Completed sampling in 1.464495 secs
```

The posterior samples of the linkage structure can be accessed by
calling `extract(fit, "links")`. However, these samples only cover the
records that were considered as candidate pairs. We can obtain samples
of the complete linkage structure (for all records) using the following
function.

``` r
links_samples <- complete_links_samples(fit, RLdata500$rec_id)
```

Since RLdata500 comes with ground truth, we can evaluate the predicted
linkage structure using supervised metrics. We find that the model
performs well on this data set, achieving high pairwise precision and
recall.

``` r
n_records <- nrow(RLdata500)
true_pairs <- membership_to_pairs(identity.RLdata500)
metrics <- apply(links_samples, 1, function(links) {
  pred_pairs <- membership_to_pairs(links)
  pr <- precision_pairs(true_pairs, pred_pairs)
  re <- recall_pairs(true_pairs, pred_pairs)
  c(Precision = pr, Recall = re)
}) %>% t()
summary(metrics)
#>    Precision          Recall      
#>  Min.   :0.8772   Min.   :0.9600  
#>  1st Qu.:0.9259   1st Qu.:1.0000  
#>  Median :0.9434   Median :1.0000  
#>  Mean   :0.9386   Mean   :0.9968  
#>  3rd Qu.:0.9615   3rd Qu.:1.0000  
#>  Max.   :0.9804   Max.   :1.0000
```

## Licence

GPL-3

## References

<div id="refs" class="references hanging-indent">

<div id="ref-fellegi1969">

Fellegi, Ivan P., and Alan B. Sunter. 1969. “A Theory for Record
Linkage.” *Journal of the American Statistical Association* 64 (328):
1183–1210. <https://doi.org/10.1080/01621459.1969.10501049>.

</div>

<div id="ref-sadinle2014">

Sadinle, Mauricio. 2014. “Detecting Duplicates in a Homicide Registry
Using a Bayesian Partitioning Approach.” *Ann. Appl. Stat.* 8 (4):
2404–34. <https://doi.org/10.1214/14-AOAS779>.

</div>

</div>
