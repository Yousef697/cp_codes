#include <bits/stdc++.h>
#define int long long

using namespace std;
using ll = long long;
using ld = long double;

struct MaxFlow {
    const int inf = 1e9 + 9;
    int n;
    vector<int> cur_path, reach;
    vector<vector<int> > cap, init;
    vector<vector<int> > adj, g;
    vector<vector<int> > paths, vis;
    vector<pair<int, int> > cuts;

    MaxFlow(int n): n(n) {
        adj = g = vector<vector<int> >(n + 1);
        cap = init = vector<vector<int> >(n + 1, vector<int>(n + 1, 0));
    }

    void add_edge(int undir, int u, int v, int c = 1) {
        init[u][v] += c;
        if (undir) init[v][u] += c;
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }
    int bfs(int s, int t, vector<int> &par) {
        par = vector<int>(n + 1, -1);
        queue<pair<int, int> > q;

        par[s] = s;
        q.emplace(s, inf);

        while (!q.empty()) {
            auto [u, flow] = q.front();
            q.pop();
            for (auto &v: adj[u]) {
                if (par[v] == -1 && cap[u][v] > 0) {
                    par[v] = u;
                    int new_flow = min(flow, cap[u][v]);
                    if (v == t) return new_flow;
                    q.emplace(v, new_flow);
                }
            }
        }
        return 0;
    }
    int max_flow(int s, int t) {
        cap = init;
        int flow = 0;
        vector<int> par;
        while (true) {
            int add = bfs(s, t, par);
            if (add == 0) break;
            flow += add;
            int cur = t;
            while (cur != s) {
                int prv = par[cur];
                cap[prv][cur] -= add;
                cap[cur][prv] += add;
                cur = prv;
            }
        }
        return flow;
    }

    void dfs_paths(int u, int p, int t) {
        vis[u][p] = 1;
        if (u == t) {
            cur_path.emplace_back(u);
            paths.emplace_back(cur_path);
            cur_path.pop_back();
            return;
        }
        cur_path.emplace_back(u);
        for (auto &v: g[u]) {
            if (v != p && !vis[v][u]) {
                dfs_paths(v, u, t);
                break;
            }
        }
        cur_path.pop_back();
    }
    int construct_paths(int s, int t) {
        int flow = max_flow(s, t);
        vis = vector<vector<int> >(n + 1, vector<int>(n + 1, 0));
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (cap[i][j] == 0 && init[i][j] == 1) {
                    g[i].emplace_back(j);
                }
            }
        }
        paths.clear();
        cur_path = {s};
        for (auto &v: g[s]) {
            dfs_paths(v, 1, t);
        }
        return flow;
    }

    void dfs_cuts(int u, int p) {
        reach[u] = 1;
        for (auto &v: adj[u]) {
            if (v != p && !reach[v] && cap[u][v])
                dfs_cuts(v, u);
        }
    }
    int construct_min_cuts(int s, int t) {
        int flow = max_flow(s, t);
        reach = vector<int>(n + 1);
        cuts.clear();
        dfs_cuts(s, s);
        for (int i = 1; i <= n; i++) {
            for (auto &j: adj[i]) {
                if (reach[i] && !reach[j])
                    cuts.emplace_back(i, j);
            }
        }
        return flow;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    return 0;
}
