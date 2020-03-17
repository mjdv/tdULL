#include "centrality.hpp"
#include <queue>
#include <stack>
#include <cmath>

std::vector<int> DegreeCentrality(const Graph &G) {
    std::vector<int> degree(G.N);
    for(int i = 0; i < G.N; i++)
        degree[i] = G.Adj(i).size();
    return degree;
}

std::vector<double> BetweennessCentrality(const Graph &G) {
    int N = G.N;

    std::vector<double> betweenness(N, 0.0);

    for(int source = 0; source < N; source++) {
        std::vector<double> delta(N, 0.0);
        std::vector<double> sigma(N, 0);
        sigma[source] = 1;

        std::vector<int> dist(N, -1);
        dist[source] = 0;
        std::vector<std::vector<int>> pred(N, std::vector<int>());
        std::queue<int> q;
        q.push(source);

        std::stack<int> s;

        while(!q.empty()) {
            int cur = q.front(); q.pop();
            s.push(cur);

            for(int nb : G.Adj(cur)) {
                if(dist[nb] == -1) {
                    dist[nb] = dist[cur] + 1;
                    q.push(nb);
                }
                if(dist[nb] == dist[cur] + 1) {
                    sigma[nb] += sigma[cur];
                    pred[nb].push_back(cur);
                }
            }
        }

        while(!s.empty()) {
            int cur = s.top(); s.pop();
            for(int back : pred[cur]) {
                delta[back] = delta[back] + (sigma[back]/sigma[cur])*(1 + delta[cur]);
            }
            if(cur != source) {
                betweenness[cur] += delta[cur];
            }
        }
    }
    return betweenness;
}

std::vector<double> EigenvectorCentrality(const Graph &G, size_t steps) {
  int N = G.N;
  std::vector<double> centrality(N, sqrt(1.0/N));
  for (size_t step = 0; step < steps; step++) {
    std::vector<double> new_centrality(N, 0);
    double norm = 0.0;
    for (size_t n = 0; n < N; n++) {
      for (auto nbr : G.Adj(n))
        new_centrality[n] += centrality[nbr];
      norm += new_centrality[n] * new_centrality[n];
    }
    norm = sqrt(norm);
    for (size_t n = 0; n < N; n++)
      new_centrality[n] /= norm;
    centrality = std::move(new_centrality);
  }
  return centrality;
}

std::vector<double> PageRankCentrality(const Graph &G, size_t steps, double damping) {
  int N = G.N;
  std::vector<double> centrality(N, sqrt(1.0/N));
  for (size_t step = 0; step < steps; step++) {
    std::vector<double> new_centrality(N, 0);
    double centrality_sum = 0.0;
    for (size_t n = 0; n < N; n++)
      centrality_sum += centrality[n];

    double norm = 0.0;
    for (size_t n = 0; n < N; n++) {
      for (auto nbr : G.Adj(n))
        new_centrality[n] += damping * centrality[nbr];
      new_centrality[n] += (1 - damping) * centrality_sum;
      norm += new_centrality[n] * new_centrality[n];
    }
    norm = sqrt(norm);
    for (size_t n = 0; n < N; n++)
      new_centrality[n] /= norm;
    centrality = std::move(new_centrality);
  }
  return centrality;
}
