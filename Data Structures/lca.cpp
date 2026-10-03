#include <bits/stdc++.h>

using namespace std;

// Lowest Common Ancestor (LCA)
struct LCA {
    int n, root, lg, cnt;
    vector<int> open, close, level;
    vector<vector<int> > table, adj;

    LCA(int n, int root = 1): n(n), root(root) {
        lg = log2(n), cnt = 0;
        open = close = level = vector<int>(n + 1, 0);
        adj = vector<vector<int>> (n + 1);
        table = vector<vector<int> >(n + 1, vector<int>(lg + 1, root));
    }
    void add_edge(int u, int v) {
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }
    void dfs(int u, int p) {
        open[u] = ++cnt;
        table[u][0] = p;
        for (int i = 1; i <= lg; i++)
            table[u][i] = table[table[u][i - 1]][i - 1];
        for (auto& v: adj[u])
            if (v != p)
                level[v] = level[u] + 1, dfs(v, u);
        close[u] = cnt;
    }
    bool is_ancestor(int u, int v) {
        return open[u] <= open[v] && close[v] <= close[u];
    }
    int query(int u, int v) {
        if (is_ancestor(u, v)) return u;
        if (is_ancestor(v, u)) return v;
        for (int i = lg; i >= 0; i--)
            if (!is_ancestor(table[u][i], v))
                u = table[u][i];
        return table[u][0];
    }
    int get_distance(int u, int v) {
        return level[u] + level[v] - 2 * level[query(u, v)];
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, q;
    cin >> n >> q;

    LCA lca(n, 1);
    for (int i = 1; i <= n - 1; i++) {
        int u, v;
        cin >> u >> v;
        lca.add_edge(u, v);
    }

    while (q--) {
        int u, v;
        cin >> u >> v;
        
        cout << lca.query(u, v) << "\n";
    }

    return 0;
}
