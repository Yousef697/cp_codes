#include <bits/stdc++.h>

using namespace std;
using ll = long long;

struct dfs_tree {
    int n, is_connected;
    vector<int> in, out, dp, vis;
    vector<vector<int> > adj, tree;
    vector<pair<int, int> > bridges, back_edges, edges;

    dfs_tree(int _n) {
        n = _n, is_connected = 1;
        in = out = dp = vis = vector<int>(n + 1, 0);
        tree = adj = vector<vector<int> >(n + 1);
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
        for (auto v: adj[u]) {
            if (v == p)
                continue;

            if (vis[v] == 1) {
                in[v]++, out[u]++;
                back_edges.push_back({u, v});
            } else if (!vis[v]) {
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

        for (auto v: tree[u]) {
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
struct UnionFind {
    int forests, N = 2e5;
    vector<int> parent, size, rank;

    UnionFind() {
        forests = N;
        parent = size = rank = vector<int>(N + 1, 0);
        for (int i = 0; i <= N; i++)
            parent[i] = i, size[i] = rank[i] = 1;
    }

    UnionFind(int n) {
        forests = n;
        parent = size = rank = vector<int>(n + 1, 0);
        for (int i = 0; i <= n; i++)
            parent[i] = i, size[i] = rank[i] = 1;
    }

    int find_set(int v) {
        if (v == parent[v])
            return v;
        return parent[v] = find_set(parent[v]);
    }

    void union_sets(int a, int b) {
        a = find_set(a);
        b = find_set(b);

        if (a == b)
            return;

        if (size[a] > size[b])
            swap(a, b);

        forests--;
        parent[a] = b;
        size[b] += size[a];
    }

    int size_set(int u) { return size[find_set(u)]; }
    bool same_set(int u, int v) { return find_set(u) == find_set(v); }
};

const int mod = 1e9 + 7;
int fp(int i, int j) {
    if (j == 0)
        return 1;
    int ret = fp(i, j /  2);
    ret = (1ll * ret * ret) % mod;
    if (j & 1)
        ret = (1ll * ret * i) % mod;
    return ret;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);


    int n, m;
    cin >> n >> m;

    dfs_tree tree(n);
    vector<pair<int, int> > e;
    for (int i = 0; i < m; i++) {
        int u, v;
        cin >> u >> v;
        tree.add_edge(u, v);
        e.push_back({u, v});
    }
    tree.work();

    vector<int> par(n + 1);
    function<void(int, int)> get_par = [&](int u, int p) -> void {
        par[u] = p;
        for (auto v: tree.tree[u]) {
            if (v == p)
                continue;
            get_par(v, u);
        }
    };
    get_par(1, 0);

    UnionFind uf(n);
    vector<int> prf(n + 1), cur(n + 1);
    for (auto &[u, v]: tree.back_edges) {
        int x = u;
        uf.union_sets(u, v);
        while (x != v) {
            x = par[x];
            uf.union_sets(x, v);
        }
        x = uf.find_set(v);
        prf[x] = 1;
        cur[x] = 1;
    }
    for (int i = 1; i <= n; i++) {
        int x = uf.find_set(i);
        // do nothing
    }

    vector<vector<int> > adj(n + 1);
    for (auto &[i, j]: e) {
        i = uf.find_set(i);
        j = uf.find_set(j);
        if (i != j) {
            adj[i].push_back(j);
            adj[j].push_back(i);
        }
    }

    int node = -1;
    for (int i = 1; i <= n; i++) {
        if (!adj[i].empty())
            node = i;
    }
    par = vector<int>(n + 1, node);
    function<void(int, int)> get_prefix = [&](int u, int p) {
        prf[u] += prf[p], par[u] = p;
        for (auto v: adj[u]) {
            if (v != p)
                get_prefix(v, u);
        }
    };
    if (~node)
        get_prefix(node, node);

    int lg = 20, time = 0;
    vector<int> in(n + 1), out(n + 1);
    vector<vector<int>> up(n + 1, vector<int>(lg, node));
    function<void(int, int)> get_lca = [&](int u, int p) {
        in[u] = ++time;
        up[u][0] = p;
        for (int i = 1; i < lg; i++)
            up[u][i] = up[up[u][i - 1]][i - 1];
        for (auto v : adj[u]) {
            if (v == p)
                continue;
            get_lca(v, u);
        }
        out[u] = time;
    };
    if (~node)
        get_lca(node, node);

    // for (int i = 1; i <= n; i++) {
    //     cout << i << ": ";
    //     for (int j = 0; j < 5; j++) {
    //         cout << up[i][j] << " ";
    //     }
    //     cout << "\n";
    // }

    auto isa = [&](int u, int v) -> bool {
        return in[u] <= in[v] && out[v] <= out[u];
    };
    auto lca = [&](int u, int v) -> int {
        if (isa(u, v)) return u;
        if (isa(v, u)) return v;

        int ans = u;
        for (int i = lg - 1; i >= 0; i--) {
            if (!isa(up[ans][i], v))
                ans = up[ans][i];
        }
        return par[ans];
    };

    // cout << node << "\n";
    // for (int i = 1; i <= n; i++) {
    //     for (auto j : adj[i]) {
    //         cout << 100 + i << " " << 100 + j << "\n";
    //     }
    // }
    // for (int i = 1; i <= n; i++) {
    //     cout << i << " " << cur[i] << "\n";
    // }
    // cout << "\n\n";
    // return 0;

    int q;
    cin >> q;

    while (q--) {
        int u, v;
        cin >> u >> v;

        // int f = u == 14 && v == 38;

        u = uf.find_set(u);
        v = uf.find_set(v);

        // cout << u << " " << v << "\n";
        // cout << 100 + u << " " << 100 + v << "\n";

        if (u == v) {
            cout << cur[u] + 1 << "\n";
            continue;
        }

        int lc = lca(u, v), p = prf[u] + prf[v] - 2 * prf[lc] + cur[lc];
        // cout << lc << " ---- ";
        if (lc == u) {
            p = prf[v] - prf[u] + cur[u];
        }
        else if (lc == v) {
            p = prf[u] - prf[v] + cur[v];
        }
        // if (f)cout << lc << "\n";
        cout << fp(2, p) << "\n";
        // cout << endl;
    }

    return 0;
}