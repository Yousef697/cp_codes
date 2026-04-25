#include <bits/stdc++.h>

using namespace std;
using ll = long long;
const int N = 2e5 + 5;

struct centroid_decomposition {
    // centroid variables
    int n;
    vector<int> sz;
    vector<bool> blocked;
    vector<vector<int>> adj;

    // problem-specific variables
    int k;
    ll ans;
    vector<int> cur;

    centroid_decomposition(int n, int k): n(n), k(k) {
        sz.assign(n + 1, 0);
        blocked.assign(n + 1, false);
        adj.assign(n + 1, {});

        ans = 0;
        cur.assign(n + 1, 0);
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

    void get_dist(int u, int p, int d, vector<int>& dist) {
        dist.push_back(d);
        if (d == k) return;
        for (auto& v : adj[u]) {
            if (v == p || blocked[v]) continue;
            get_dist(v, u, d + 1, dist);
        }
    }

    void process_centroid(int centroid) {

        cur[0] = 1;
        vector<int> all_dist;
        for (auto& v : adj[centroid]) {
            if (blocked[v]) continue;

            vector<int> dist;
            get_dist(v, centroid, 1, dist);
            for (auto& d : dist) {
                ans += cur[k - d];
            }
            for (auto& d : dist) {
                cur[d]++;
                all_dist.push_back(d);
            }
        }
        for (auto& d : all_dist) {
            cur[d] = 0;
        }
        cur[0] = 0;
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
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, k;
    cin >> n >> k;

    centroid_decomposition cd(n, k);

    for (int i = 1; i < n; i++) {
        int u, v;
        cin >> u >> v;
        cd.add_edge(u, v);
    }

    cd.solve();
    cout << cd.ans << "\n";

    return 0;
}