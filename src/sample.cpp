#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <queue>
#include "updates.h"

Rcpp::NumericVector unlist_vec(const Rcpp::List& list)
{
  std::size_t out_size = 0;
  for (const auto& x : list) { out_size += Rf_length(x); }
  Rcpp::NumericVector out(out_size);

  std::size_t i = 0;
  for (const Rcpp::NumericVector& x : list)
  {
    std::copy(x.begin(), x.end(), out.begin() + i);
    i += x.size();
  }

  return out;
}

//[[Rcpp::export(.sample)]]
Rcpp::S4 sample(const Rcpp::S4& init_state, int n_samples, int thin_interval=1, 
                int burnin_interval=0) 
{
  int n_iter = burnin_interval + thin_interval * n_samples;
  Progress p(n_iter, true);
  
  // Don't want to overwrite input state
  Rcpp::S4 state = Rcpp::clone(init_state);
  
  // Preallocate R arrays to store history of variables along the chain
  const Rcpp::List& priors = state.slot("priors");
  const Rcpp::List& lambda = priors["lambda"];
  std::size_t n_attributes = lambda.size();
  Rcpp::IntegerVector n_levels = Rcpp::sapply(lambda, Rf_length);
  int ncol = Rcpp::sum(n_levels);
  const Rcpp::IntegerVector& rec_ids = state.slot("rec_ids");
  Rcpp::IntegerVector hist_n_clusters(n_samples);
  Rcpp::IntegerMatrix hist_links(n_samples, rec_ids.size());
  Rcpp::NumericMatrix hist_m(n_samples, ncol);
  Rcpp::NumericMatrix hist_u(n_samples, ncol);
  colnames(hist_links) = rec_ids;
  
  const Rcpp::CharacterVector& attr_names = lambda.names();
  Rcpp::CharacterVector m_colnames(ncol);
  std::size_t a, l, i=0;
  for (a = 0; a < n_attributes; ++a) {
    for (l = 0; l < Rf_length(lambda[a]); ++l) {
      m_colnames[i] = Rcpp::as<std::string>(attr_names[a]) + "[" + std::to_string(l + 1) + "]";
      ++i;
    }
  }
  colnames(hist_m) = m_colnames;
  colnames(hist_u) = m_colnames;

  int sample_ctr = 0;
  int completed_iter = Rcpp::as<Rcpp::IntegerVector>(state.slot("iteration"))[0];  
  // number of completed iterations (may be different from state's iteration counter if resuming)
  
  // Build inverted index (from clust_id to rec_ids assigned to the cluster)
  const Rcpp::IntegerVector& links = state.slot("links");
  Map<int, Set<int>> links_inv;
  int rec_id = 0;
  for (const int& clust_id : links) {
    links_inv[clust_id].insert(rec_id);
    rec_id++;
  }
  
  // Queue to allow recycling of empty cluster ids
  std::queue<int> empty_clust_ids;
  int clust_id = 0;
  for (; clust_id < links.size(); clust_id++) {
    auto it = links_inv.find(clust_id);
    if (it == links_inv.end()) {
      // No entry for this clust_id, so must be empty
      empty_clust_ids.push(clust_id);
    }
  }
  
  while (sample_ctr < n_samples) 
  {
    // Update unobserved variables
    update_m_u_probs(state);
    update_links(state, links_inv, empty_clust_ids);
    
    completed_iter++;
    state.slot("iteration") = completed_iter;
    
    if (completed_iter >= burnin_interval) 
    {
      // Finished burn-in, so start saving samples
      if ((completed_iter - burnin_interval) % thin_interval == 0)
      {
        // Update history using this sample
        hist_n_clusters[sample_ctr] = links_inv.size();
        const Rcpp::List& links = state.slot("links");
        hist_links.row(sample_ctr) = links;
        const Rcpp::List& m = state.slot("m");
        const Rcpp::List& u = state.slot("u");
        Rcpp::NumericVector m_vec = unlist_vec(m);
        Rcpp::NumericVector u_vec = unlist_vec(u);
        hist_m.row(sample_ctr) = m_vec;
        hist_u.row(sample_ctr) = u_vec;
        sample_ctr++;
      }
    }
    
    if (Progress::check_abort()) { break; }
    p.increment();
  }
  
  // Prepare output
  Rcpp::S4 result("BDDFit");
  Rcpp::List history = Rcpp::List::create(
    Rcpp::Named("n_clusters") = hist_n_clusters,
    Rcpp::Named("links") = hist_links,
    Rcpp::Named("m") = hist_m,
    Rcpp::Named("u") = hist_u
  );
  result.slot("history") = history;
  result.slot("state") = state;
  return result;
}

