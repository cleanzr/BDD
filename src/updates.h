#pragma once

#include <Rcpp.h>
#include <unordered_map>
#include <set>
#include <queue>

template <class K, class V>
using Map = std::unordered_map<K, V>;

template <class V>
using Set = std::set<V>;

template <class V>
using Vec = std::vector<V>;

void update_m_u_probs(Rcpp::S4& state);
void update_links(Rcpp::S4& state, Map<int, Set<int>>& links_inv, 
                  std::queue<int>& empty_clust_ids);