#include <bits/stdc++.h>

using namespace std;

// Lowest Common Ancestor (LCA)
struct LCA {
    /// @brief LCA Variable
    int lg = log2(2e5 + 5), cnt = 1, n;
    vector<int> open, close, level;
    vector<vector<int> > table, adj;

    /// @brief Contructor
    LCA(int root, vector<vector<int> > &a) {
        n = a.size() - 1;
        lg = log2(n);
        open = close = level = vector<int>(n + 1, 0);
        table = vector<vector<int> >(n + 1, vector<int>(lg + 1, root));
        adj = a;

        dfs(root, root);
    }

    /// @brief do depth-first search on the tree to know if a node is an ancestor of another node
    /// @param node
    /// @param parent
    void dfs(int node, int parent) {
        open[node] = cnt++;

        table[node][0] = parent;
        for (int i = 1; i <= lg; i++) {
            int x = table[node][i - 1];
            table[node][i] = table[x][i - 1];
        }

        for (auto i: adj[node]) {
            if (i != parent)
                level[i] = level[node] + 1, dfs(i, node);
        }

        close[node] = cnt++;
    }

    /// @brief check if a node "a" is an ancestor of node "b"
    /// @param a
    /// @param b
    /// @return true if node "a" is an ancestor of node "b", false otherwise
    bool is_ancestor(int a, int b) {
        return open[a] <= open[b] && close[b] <= close[a];
    }

    /// @brief get the least common ancestor of two nodes
    /// @param a
    /// @param b
    /// @return LCA
    int query(int a, int b) {
        if (is_ancestor(a, b))
            return a;
        if (is_ancestor(b, a))
            return b;

        for (int i = lg; i >= 0; i--) {
            if (!is_ancestor(table[a][i], b)) {
                a = table[a][i];
            }
        }

        return table[a][0];
    }

    /// @brief get the distance between two nodes in a tree
    /// @param a
    /// @param b
    /// @return distance between node "a", and "b"
    int get_distance(int u, int v) {
        int ancestor = query(u, v);
        return abs(level[ancestor] - level[v]) + abs(level[ancestor] - level[u]);
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, q;
    cin >> n >> q;

    vector<vector<int> > adj(n + 1);
    for (int i = 1; i <= n - 1; i++) {
        int u, v;
        cin >> u >> v;

        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    LCA lca(1, adj);

    while (q--) {
        int u, v;
        cin >> u >> v;

        cout << lca.query(u, v) << "\n";
    }

    return 0;
}
