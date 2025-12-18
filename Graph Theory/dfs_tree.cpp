#include <bits/stdc++.h>

using namespace std;
using ll = long long;

struct dfs_tree {
    int n, is_connected;
    vector<int> in, out, dp, vis;
    vector<vector<int>> adj, tree;
    vector<pair<int, int>> bridges, back_edges, edges;

    dfs_tree(int _n) {
        n = _n, is_connected = 1;
        in = out = dp = vis = vector<int>(n + 1, 0);
        tree = adj = vector<vector<int>>(n + 1);
    }

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }
    void work() {
        dfs(1, 0);
        for (int i = 1; i <= n; i++)
            if (!vis[i]) is_connected = 0;
        vis = vector<int>(n + 1);
        solve(1, 0);
    }
    void dfs(int u, int p) {
        vis[u] = 1;
        for (auto v : adj[u]) {
            if (v == p)
                continue;

            if (vis[v] == 1) {
                in[v]++, out[u]++;
                back_edges.push_back({u, v});
            }
            else if (!vis[v]) {
                tree[u].push_back(v);
                tree[v].push_back(u);
                edges.push_back({u, v});
                dfs(v, u);
            }
        }
        vis[u] = 2;
    }
    void solve(int u, int p) {
        vis[u] = 1;
        int ok = out[u] - in[u];

        for (auto v : tree[u]) {
            if (v == p)
                continue;
            solve(v, u);
            ok += dp[v];
        }

        dp[u] = ok;
        if (ok == 0) {
            if (u != 1 || adj[u].size() == 1) {
                bridges.push_back({p, u});
            }
        }
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    int n, m;
    cin >> n >> m;

    dfs_tree tree(n);
    for (int i = 1; i <= m; i++) {
        int u, v;
        cin >> u >> v;
        tree.add_edge(u, v);
    }

    tree.work();

    if (!tree.bridges.empty() || !tree.is_connected) {
        cout << "IMPOSSIBLE\n";
        return 0;
    }

    for (auto& [u, v] : tree.edges)
        cout << u << " " << v << "\n";
    for (auto& [u, v] : tree.back_edges)
        cout << u << " " << v << "\n";

    return 0;
}
