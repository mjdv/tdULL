#include "centrality.hpp"
#include <queue>
#include <stack>

std::vector<int> DegreeCentrality(const SubGraph &G) {
    std::vector<int> degree(G.vertices.size());
    for(int i = 0; i < G.vertices.size(); i++)
        degree[i] = G.Adj(i).size();
    return degree;
}

std::vector<double> BetweennessCentrality(const SubGraph &G) {
    int N = G.vertices.size();

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
