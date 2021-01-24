#include <Rcpp.h>
#include <algorithm>
#include "updates.h"

//' Helper function for `get_clust_weights`
Vec<Vec<double> > compute_log_m_on_u(const Rcpp::S4& state) 
{
  const Rcpp::List& m = state.slot("m");
  const Rcpp::List& u = state.slot("u");
  
  std::size_t n_attributes = m.length();
  std::size_t n_levels, a, i;

  Vec<Vec<double>> out(n_attributes);
  for (a = 0; a < n_attributes; ++a) {
    const Rcpp::NumericVector& this_m = m[a];
    const Rcpp::NumericVector& this_u = u[a];
    
    n_levels = this_m.size() + 1;
    Vec<double> this_log_m_on_u(n_levels);
    double log_m_1, log_m_2 = 0.0, log_u_1, log_u_2 = 0.0;
    for (i = 0; i < n_levels; ++i) {
      log_m_1 = (i == n_levels - 1) ? 0.0 : std::log(this_m[i]);
      log_m_2 += (i == 0) ? 0.0 : std::log(1 - this_m[i-1]);
      log_u_1 = (i == n_levels - 1) ? 0.0 : std::log(this_u[i]);
      log_u_2 += (i == 0) ? 0.0 : std::log(1 - this_u[i-1]);
      this_log_m_on_u[i] = log_m_1 + log_m_2 - log_u_1 - log_u_2;
    }
    out[a] = this_log_m_on_u;
  }
  return out;
}

//' Helper function for `update_links`
//'
//' Computes the pmf over valid cluster ids for a particular record id
Vec<double> get_clust_weights(int rec_id, const Set<int> &valid_clust_ids, 
                              const Map<int, Set<int> > &links_inv, 
                              const Rcpp::S4& state)
{
  const Rcpp::List& row_id_index = state.slot("row_id_index");
  const Rcpp::List& comparisons = state.slot("comparisons");
  Vec<Vec<double> > log_m_on_u = compute_log_m_on_u(state);
  
  // Extra entry at the end corresponds to an empty cluster. It will be assigned a 
  // weight of 1.0.
  Vec<double> weights(valid_clust_ids.size() + 1, 0.0);
  const Rcpp::IntegerVector& row_ids = row_id_index[rec_id];

  std::size_t i = 0, a;
  std::size_t n_attributes = comparisons.size();

  double lambda;
  int comp_val;
  for (const auto &clust_id : valid_clust_ids) {
    auto it = links_inv.find(clust_id);
    if (it != links_inv.end()) {
      const Set<int>& cluster = it->second;
      for (const auto& rec_id_pr : cluster) {
        const Rcpp::IntegerVector& row_ids_pr = row_id_index[rec_id_pr];
        Vec<int> result;
        std::set_intersection(row_ids.begin(), row_ids.end(), 
                              row_ids_pr.begin(), row_ids_pr.end(), 
                              std::back_inserter(result));
        for (const auto& row_id : result) {
          // Compute lambda
          lambda = 0.0;
          // Loop over attributes
          for (a = 0; a < n_attributes; ++a) {
            // If observed
            const Rcpp::IntegerVector& comp_vals = comparisons[a];
            comp_val = comp_vals[row_id]; 
            if (comp_val != NA_INTEGER) {
              lambda += log_m_on_u[a][comp_val - 1]; // R indexing
            }
          }
          //Rcpp::Rcout << "row_id=" << row_id << " lambda=" << lambda << "\n";
          weights[i] += lambda;
        }
      }
    }
    i++;
  }
  
  // Exponentiate to convert from log-likelihood to likelihood
  for (auto& w : weights) { w = exp(w); }

  return weights;
}

//' Helper function for `update_links`
//'
//' Computes valid cluster ids given a list of possible pairing record ids
Set<int> get_valid_clust_ids(const Rcpp::IntegerVector& pssbl_rec_ids, 
                             const Rcpp::IntegerVector& links, 
                             const Map<int, Set<int> >& links_inv)
{
  bool valid_clust = true;

  // Start with the set of clusters currently occupied by at least one of the 
  // possible pairing records.
  Set<int> pssbl_clust_ids;
  for (const auto& rec_id : pssbl_rec_ids) pssbl_clust_ids.insert(links[rec_id]);

  // Iterate over these clusters, removing any that are invalid
  for (auto it = pssbl_clust_ids.begin(); it != pssbl_clust_ids.end(); ) {
    // Need to look at rec_ids in the cluster
    auto links_it = links_inv.find(*it);
    if (links_it != links_inv.end()) {
      const Set<int>& cluster = links_it->second;
      // Cluster is valid if all the rec_ids appear in pssbl_rec_ids
      Set<int> diff;
      std::set_difference(cluster.begin(), cluster.end(), 
                          pssbl_rec_ids.begin(), pssbl_rec_ids.end(), 
                          std::inserter(diff, diff.end()));
      if (!diff.empty()) valid_clust = false;
    }
    
    // Increment iterator correctly when deleting an element
    if (!valid_clust) {
      pssbl_clust_ids.erase(it++);
      valid_clust = true;
    } else {
      ++it;
    }
  }

  return pssbl_clust_ids;
}

//' Update the Links
//' 
//' @param state[in,out] A BDDModel class representing the current state of 
//'   the Markov chain.
//' @param links_inv[in,out] Inverted index for the linkage structure, where 
//'   the keys are the cluster ids and the values are the record ids currently 
//'   assigned to that cluster.
//' @param empty_clust_ids[in,out] Queue which keeps track of cluster ids that 
//'   are unused (have no records assigned to them).
void update_links(Rcpp::S4& state, Map<int, Set<int>>& links_inv, 
                  std::queue<int>& empty_clust_ids)
{ 
  int rec_id, old_clust_id, new_clust_id, n_valid_clust_ids;
  const Rcpp::IntegerVector& links_old = state.slot("links");
  const Rcpp::List& pair_index = state.slot("pair_index");
  Rcpp::IntegerVector links = Rcpp::clone(links_old);
  
  for(rec_id = 0; rec_id < links.size(); rec_id++){
    // Remove from current cluster
    old_clust_id = links[rec_id];
    links[rec_id] = -1;
    auto it = links_inv.find(old_clust_id);
    if (it != links_inv.end()) { 
      it->second.erase(rec_id); 
      if (it->second.empty()) {
        // Add to queue of empty clusters
        empty_clust_ids.push(old_clust_id);
        // Remove from inverted index
        links_inv.erase(it);
      }
    }
    
    const Rcpp::IntegerVector& pssbl_rec_ids = pair_index[rec_id];
    
    /* assumes that pssbl_rec_ids is sorted */
    // Filters out any record ids in pssbl_rec_ids that are invalid
    Set<int> valid_clust_ids = get_valid_clust_ids(pssbl_rec_ids, links, links_inv);

    if (valid_clust_ids.empty()) {
      new_clust_id = empty_clust_ids.front();
    } else {
      // Get weight for each valid cluster (last entry corresponds to a "new" 
      // empty cluster)
      Vec<double> weights = get_clust_weights(rec_id, valid_clust_ids, links_inv, state);
      
      // Sample
      n_valid_clust_ids = valid_clust_ids.size();
      Rcpp::NumericVector weights_R(weights.begin(), weights.end());
      int set_idx = Rcpp::sample(n_valid_clust_ids + 1, 1, true, weights_R, false)[0];
      if (set_idx < n_valid_clust_ids) {
        // An existing cluster (need to do sequential access for set)
        new_clust_id = *std::next(valid_clust_ids.begin(), set_idx);
      } else {
        // A "new" empty cluster
        new_clust_id = empty_clust_ids.front();
      }
    }
    
    // Assign record to new cluster
    links[rec_id] = new_clust_id;
    
    // Try to create entry in inverted index assuming cluster doesn't already exist
    auto result = links_inv.insert(std::make_pair(new_clust_id, Set<int> {rec_id}));
    
    // If we didn't insert anything, cluster already exists in the index. 
    // Need to insert in existing set.
    if (!result.second) { (result.first)->second.insert(rec_id); }
    
    // Update empty_clust_ids
    if (new_clust_id == empty_clust_ids.front()) { empty_clust_ids.pop(); }
  }
  
  state.slot("links") = links;
}

//' Tabulation for an Integer Vector
//' 
//' Computes `tabulate(x, nbins)` (R notation) for an integer vector `x`.
//' 
//' @param x Data vector with integer values in the range [1, nbins]
//' @param nbins Specifies the range of the integer values
//' @return An integer vector, where the i-th entry counts the number of 
//'   values in x equal to i (assuming 1-based indexing)
Vec<int> tabulate(const Rcpp::IntegerVector& x, int nbins) 
{
  Vec<int> counts(nbins, 0);
  for (auto& xi : x) 
  {
    if (xi > 0 && xi <= nbins)
      counts[xi - 1]++;
  }
  return counts;
}

//' Specialized Cumulative Sum
//' 
//' Computes `rev(cumsum(rev(x)))[-1]` (in R notation). This quanity is 
//' required when updating the m and u probabilities.
Vec<int> compute_cumsum(const Vec<int>& x) 
{
  Vec<int> out(std::max(int(x.size() - 1), 0));
  std::partial_sum(x.rbegin(), std::prev(x.rend()), out.begin());
  std::reverse(out.begin(), out.end());
  return out;
}

//' Update m probability
//' 
//' The conditional posterior for `m` is a truncated Beta distribution. This 
//' function implements an auxiliary variable sampling scheme to update m 
//' as described by Damien & Walker (2001).
//' 
//' @param m Current value of `m`. This is needed as an auxiliary variable 
//'   `y` is introduced to update `m` by drawing y | m then m | y.
//' @param shape1,shape2 Non-negative shape parameters of the Beta distribution
//' @param a Lower truncation point. Must be in the interval [0, 1).
//' @return New value of m.
//' 
//' @references
//' P. Damien & S. G. Walker (2001). Sampling Truncated Normal, Beta, and Gamma 
//' Densities, Journal of Computational and Graphical Statistics, 10:2, 
//' 206-215, DOI: 10.1198/10618600152627906. 
double update_m(double m, double shape1, double shape2, double a) 
{
  if (shape1 <= 0 || shape2 <= 0) Rcpp::stop("shape parameters must be positive");
  if (a < 0 || a >= 1) Rcpp::stop("`a` must be on the interval [0, 1)");
  
  // Posterior is not truncated
  if (a == 0) return R::rbeta(shape1, shape2);
  
  double u1, a_to_shape1;
  double m_new;
  
  a_to_shape1 = std::pow(a, shape1);
  
  if (shape2 == 1) 
  {
    // Don't need to introduce auxiliary variable. Just use inverse CDF method.
    u1 = R::unif_rand();
    m_new = std::pow(u1 * (1 - a_to_shape1) + a_to_shape1, 1/shape1);
  } 
  else 
  {
    // Need to introduce auxiliary variable y. See p. 210 in Damien & Walker 
    // (2001).
    double u2, yp, max_a_yp_to_shape1;
    
    // Draw y | m. The conditional distribution is uniform on the interval 
    // [0, (1 - m)^(shape2 - 1)]. Rather than computing y directly, we instead 
    // compute y' = 1 - y^(1/(shape2 - 1)), since it is what appears in the 
    // conditional density for m | y below.
    u1 = R::unif_rand();
    yp = std::abs(1.0 - std::pow(u1, 1/(shape2 - 1)) * (1.0 - m));
    
    // Draw m | y using the inverse CDF method. 
    u2 = R::unif_rand();
    if (shape2 >= 1) {
      m_new = std::pow(u2 * std::pow(yp, shape1) + (1 - u2) * a_to_shape1, 1/shape1);  
    } else {
      max_a_yp_to_shape1 = std::pow(std::max(a, yp), shape1);
      m_new = std::pow(u2 + (1 - u2) * max_a_yp_to_shape1, 1/shape1);
    }
  }
  
  if (m_new < a) {
    Rcpp::Function warning("warning");
    warning("update_m: new value of m=" + std::to_string(m_new) + 
      " falls below lower truncation point a=" + std::to_string(a) + 
      ". Setting m=a." );
    m_new = a;
  }
  
  return m_new;
}

//' Update the m and u Probabilities 
//' 
//' @param state[in,out] A BDDModel class representing the current state of 
//'   the Markov chain.
void update_m_u_probs(Rcpp::S4& state) 
{
  const Rcpp::IntegerVector& links = state.slot("links");
  const Rcpp::List& c_pairs = state.slot("candidate_pairs");
  const Rcpp::List& comparisons = state.slot("comparisons");
  const Rcpp::List& level_counts = state.slot("level_counts");
  const Rcpp::IntegerVector& c_pairs_1 = c_pairs["ID_1"];
  const Rcpp::IntegerVector& c_pairs_2 = c_pairs["ID_2"];
  
  // Determine which record pairs are assigned to the same cluster
  Rcpp::LogicalVector same_clust(c_pairs_1.size());
  std::size_t n_same_clust = 0;
  std::size_t i, a, l;

  bool this_same_clust;
  std::size_t n_candidates = c_pairs_1.size();
  for (i = 0; i < n_candidates; i++) {
    this_same_clust = links[c_pairs_1[i] - 1] == links[c_pairs_2[i] - 1];
    n_same_clust += this_same_clust;
    same_clust[i] = this_same_clust;
  }
  
  std::size_t n_attributes = level_counts.size();
  Vec<Vec<int>> A1(n_attributes);
  Vec<Vec<int>> A0(n_attributes);
  for (a = 0; a < n_attributes; ++a) {
    const Rcpp::IntegerVector& this_lvl_counts = level_counts[a];
    std::size_t n_lvls = this_lvl_counts.size();
    if (n_same_clust == 0) {
      A1[a] = Vec<int>(n_lvls, 0);
      A0[a] = Rcpp::as<Vec<int>>(this_lvl_counts);
    } else {
      const Rcpp::IntegerVector& cmp_values = comparisons[a];
      Vec<int> this_A1 = tabulate(cmp_values[same_clust], n_lvls);
      Vec<int> this_A0(n_lvls);
      
      std::transform(this_lvl_counts.begin(), this_lvl_counts.end(), 
                     this_A1.begin(), this_A0.begin(), std::minus<int>());
      A1[a] = this_A1;
      A0[a] = this_A0;
    }
  }
  
  Vec<Vec<int>> cum_A1(n_attributes);
  Vec<Vec<int>> cum_A0(n_attributes);
  for (a = 0; a < n_attributes; ++a) {
    // Compute sum_{h > l} a^1_{fh}
    cum_A1[a] = compute_cumsum(A1[a]);
    // Compute sum_{h > l} a^0_{fh}
    cum_A0[a] = compute_cumsum(A0[a]);
  }
  
  // Remove entry for last level of agreement
  for (a = 0; a < n_attributes; ++a) {
    A1[a].pop_back();
    A0[a].pop_back();
  }
  
  // Draw from posterior
  double m_a, m_b, u_a, u_b;
  const Rcpp::List& priors = state.slot("priors");
  const Rcpp::List& nc_counts = state.slot("nc_counts");
  const Rcpp::List& alpha0 = priors["alpha0"];
  const Rcpp::List& beta0 = priors["beta0"];
  const Rcpp::List& alpha1 = priors["alpha1"];
  const Rcpp::List& beta1 = priors["beta1"];
  const Rcpp::List& nc_A0 = nc_counts["A0"];
  const Rcpp::List& nc_cum_A0 = nc_counts["cum_A0"];
  const Rcpp::List& lambda = priors["lambda"];
  const Rcpp::List& m_old = state.slot("m");
  
  Rcpp::List m(n_attributes);
  Rcpp::List u(n_attributes);
  for (a = 0; a < n_attributes; ++a) {
    std::size_t len_A1 = A1[a].size();
    Rcpp::NumericVector this_m(len_A1);
    Rcpp::NumericVector this_u(len_A1);
    const Rcpp::NumericVector& this_m_old = m_old[a];
    const Rcpp::NumericVector& this_alpha1 = alpha1[a];
    const Rcpp::NumericVector& this_beta1 = beta1[a];
    const Rcpp::NumericVector& this_alpha0 = alpha0[a];
    const Rcpp::NumericVector& this_beta0 = beta0[a];
    const Rcpp::NumericVector& this_nc_A0 = nc_A0[a];
    const Rcpp::NumericVector& this_nc_cum_A0 = nc_cum_A0[a];
    const Rcpp::NumericVector& this_lambda = lambda[a];
    for (l = 0; l < len_A1; l++) {
      // Compute posterior shape parameters by adding prior
      m_a = A1[a][l] + this_alpha1[l];
      m_b = cum_A1[a][l] + this_beta1[l];
      u_a = A0[a][l] + this_nc_A0[l] + this_alpha0[l];
      u_b = cum_A0[a][l] + this_nc_cum_A0[l] + this_beta0[l];
      this_m[l] = update_m(this_m_old[l], m_a, m_b, this_lambda[l]);
      this_u[l] = R::rbeta(u_a, u_b);
    }
    m[a] = this_m;
    u[a] = this_u;
  }
  
  state.slot("m") = m;
  state.slot("u") = u;
}
