#include <bits/stdc++.h>

using namespace std;

// Max Flow, Edmonds-Karp, Ford-Fulkerson
/**
    Maximum Disjoint Paths: Set all edges weights by 1
    Vertex Splitting: in(i) = i, out(i) = i + n, cost of node i is c => put edge(out(i), to(i), c)
    Maximum Independent Paths: Maximum Disjoint Paths + Vertex Splitting
    Multi-Source Multi-Sink Graph: Put a Super Source and connect it to all sources (the same to sinks)
    Maximum Bipartite Graph Matching
    Minimum Paths Coverage in DAG
    Minimum Edges Coverage
    Max Independent Set
    Minimum Cut Edges
**/
vector<vector<int>> adj;
vector<vector<long long>> adj_mat;
long long Edmonds_Karp(int st, int en)
{
    int id = 0;
    int n = adj.size() - 1;
    long long max_flow = 0;
    vector<int> vis(n + 1, -1);
    vector<int> parent(n + 1, -1);

    function<long long(void)> bfs = [&]()
    {
        queue<pair<int, long long>> q;
        q.push({st, 1e18});
        parent[st] = st;

        while (!q.empty())
        {
            auto [u, w] = q.front();
            q.pop();

            for (int v : adj[u])
            {
                if (vis[v] != id && adj_mat[u][v])
                {
                    vis[v] = id;
                    parent[v] = u;
                    q.push({v, min(1LL * w, adj_mat[u][v])});

                    if (v == en)
                        return min(1LL * w, adj_mat[u][v]);
                }
            }
        }

        return 0LL;
    };

    while (true)
    {
        id++;
        long long attemp = bfs();
        if (attemp == 0)
            break;

        max_flow += attemp;

        int x = en;
        while (x != st)
        {
            int y = parent[x];
            adj_mat[y][x] -= attemp;
            adj_mat[x][y] += attemp;
            x = y;
        }
    }

    return max_flow;
}

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}