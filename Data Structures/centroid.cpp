#include <bits/stdc++.h>

using namespace std;
using ll = long long;
const int inf = 1e9 + 5;

struct centroid_decomposition {
    // centroid variables
    int n;
    vector<int> sz;
    vector<bool> blocked;
    vector<vector<int>> adj;

    // problem-specific variables
    vector<int> ans;
    vector<vector<pair<int, int>>> up;

    centroid_decomposition(int n): n(n) {
        sz.assign(n + 1, 0);
        blocked.assign(n + 1, false);
        adj.assign(n + 1, {});

        ans.assign(n + 1, inf);
        up.assign(n + 1, {});
    }

    void add_edge(int u, int v) {
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    void pre(int u, int p) {
        sz[u] = 1;
        for (auto& v : adj[u]) {
            if (v == p || blocked[v]) continue;
            pre(v, u);
            sz[u] += sz[v];
        }
    }

    int find_centroid(int u, int p, int total) {
        for (auto& v : adj[u]) {
            if (v == p || blocked[v]) continue;
            if (2 * sz[v] > total)
                return find_centroid(v, u, total);
        }
        return u;
    }

    void get_dist(int u, int p, int d, int cent) {
        up[u].push_back({cent, d});
        for (auto& v : adj[u]) {
            if (v == p || blocked[v]) continue;
            get_dist(v, u, d + 1, cent);
        }
    }

    void process_centroid(int centroid) {
        get_dist(centroid, -1, 0, centroid);
    }

    void solve(int u = 1) {
        pre(u, -1);

        int centroid = find_centroid(u, -1, sz[u]);
        process_centroid(centroid);
        blocked[centroid] = true;

        for (auto& v : adj[centroid]) {
            if (blocked[v]) continue;
            solve(v);
        }
    }

    void update(int u) {
        for (auto& [cent, dist] : up[u]) {
            ans[cent] = min(ans[cent], dist);
        }
    }

    int query(int u) {
        int ret = ans[u];
        for (auto& [cent, dist] : up[u]) {
            ret = min(ret, ans[cent] + dist);
        }
        return ret;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    int n, m;
    cin >> n >> m;

    centroid_decomposition cd(n);
    for (int i = 1; i < n; i++) {
        int u, v;
        cin >> u >> v;
        cd.add_edge(u, v);
    }

    cd.solve();
    cd.update(1);

    while (m--) {
        int t, u;
        cin >> t >> u;

        if (t == 1) cd.update(u);
        else cout << cd.query(u) << "\n";
    }
    
    return 0;
}
