#include <bits/stdc++.h>
#define int long long

using namespace std;
using ll = long long;
using ld = long double;

// Ford-Fulkerson Edmond-Karp
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
        function<void(int, int)> dfs_paths = [&](int u, int p) {
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
                    dfs_paths(v, u);
                    break;
                }
            }
            cur_path.pop_back();
        };
        for (auto &v: g[s]) {
            dfs_paths(v, 1);
        }
        return flow;
    }
    int construct_min_cuts(int s, int t) {
        int flow = max_flow(s, t);
        reach = vector<int>(n + 1);
        cuts.clear();
        function<void(int, int)> dfs_cuts = [&](int u, int p) {
            reach[u] = 1;
            for (auto &v: adj[u]) {
                if (v != p && !reach[v] && cap[u][v])
                    dfs_cuts(v, u);
            }
        };
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

// ==========================================================================================

// Push-Relabel Algorithm
struct PushRelabel {
    const int inf = 1e15;
    int n;
    queue<int> excess_vertices;
    vector<pair<int, int>> cuts;
    vector<int> excess, height, next_son;
    vector<vector<int>> cap, flow, init, paths;

    PushRelabel(int n): n(n) {
        cap = flow = init = vector<vector<int>>(n + 1, vector<int>(n + 1));
        excess = height = next_son = vector<int>(n + 1);
    }

    void add_edge(int u, int v, int c = 1) {
        init[u][v] += c;
    }
    void relabel(int u) {
        int d = inf;
        for (int v = 1; v <= n; v++) {
            if (cap[u][v] > flow[u][v])
                d = min(d, height[v]);
        }
        if (d < inf)
            height[u] = d + 1;
    }
    void push(int u, int v) {
        int d = min(excess[u], cap[u][v] - flow[u][v]);
        flow[u][v] += d;
        flow[v][u] -= d;
        excess[u] -= d;
        excess[v] += d;
        if (d && excess[v] == d)
            excess_vertices.emplace(v);
    }
    void discharge(int u) {
        while (excess[u] > 0) {
            if (next_son[u] <= n) {
                int v = next_son[u];
                if (cap[u][v] > flow[u][v] && height[u] > height[v]) {
                    push(u, v);
                }
                else {
                    next_son[u]++;
                }
            }
            else {
                relabel(u);
                next_son[u] = 1;
            }
        }
    }
    int max_flow(int s, int t) {
        cap = init;
        height[s] = n;
        excess[s] = inf;

        for (int i = 1; i <= n; i++) {
            if (i != s)
                push(s, i);
        }
        while (!excess_vertices.empty()) {
            int u = excess_vertices.front();
            excess_vertices.pop();
            if (u != s && u != t) {
                discharge(u);
            }
        }
        int ans = 0;
        for (int i = 1; i <= n; i++) {
            ans += flow[i][t];
        }
        return ans;
    }

    int construct_paths(int s, int t) {
        int ans = max_flow(s, t);
        vector<int> cur = {s};
        vector<vector<int>> vis(n + 1, vector<int>(n + 1));
        function<void(int, int)> dfs = [&](int u, int p) {
            vis[u][p] = 1;
            if (u == t) {
                cur.emplace_back(u);
                paths.emplace_back(cur);
                cur.pop_back();
                return;
            }
            cur.emplace_back(u);
            for (int v = 1; v <= n; v++) {
                if (flow[u][v] > 0 && !vis[v][u]) {
                    dfs(v, u);
                    break;
                }
            }
            cur.pop_back();
        };
        for (int i = 1; i <= n; i++) {
            if (flow[s][i] > 0)
                dfs(i, s);
        }
        return ans;
    }
    int min_cut(int s, int t) {
        int ans = max_flow(s, t);
        vector<int> vis(n + 1);
        function<void(int, int)> dfs = [&](int u, int p) {
            vis[u] = 1;
            for (int v = 1; v <= n; v++) {
                if (!vis[v] && cap[u][v] > flow[u][v])
                    dfs(v, u);
            }
        };
        dfs(s, 0);
        for (int i = 1; i <= n; i++) {
            for (int j = 1; j <= n; j++) {
                if (init[i][j] && vis[i] && !vis[j]) {
                    cuts.emplace_back(i, j);
                }
            }
        }
        return ans;
    }
};

// ==========================================================================================

// Hopcroft-Karp
struct HopcroftKarp {
    const int inf = 1e9 + 9;
    int n, m;
    vector<int> match, dist;
    vector<vector<int>> adj;
    vector<pair<int, int>> edges;

    HopcroftKarp(int n, int m): n(n), m(m) {
        match = dist = vector<int>(n + m + 1);
        adj = vector<vector<int>>(n + m + 1);
    }
    void add_edge(int u, int v) {
        v += n;
        adj[u].emplace_back(v);
        adj[v].emplace_back(u);
    }
    bool bfs() {
        queue<int> q;
        for (int i = 1; i <= n; i++) {
            if (!match[i])
                dist[i] = 0, q.emplace(i);
            else
                dist[i] = inf;
        }
        dist[0] = inf;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            if (dist[u] >= dist[0]) continue;
            for (auto& v : adj[u]) {
                if (dist[match[v]] == inf) {
                    dist[match[v]] = dist[u] + 1;
                    q.emplace(match[v]);
                }
            }
        }
        return dist[0] != inf;
    }
    bool dfs(int u) {
        if (u == 0) return true;
        for (auto& v : adj[u]) {
            if (dist[match[v]] == dist[u] + 1 && dfs(match[v])) {
                match[u] = v;
                match[v] = u;
                return true;
            }
        }
        dist[u] = inf;
        return false;
    }
    int calc() {
        int ans = 0;
        while (bfs()) {
            for (int i = 1; i <= n; i++) {
                if (!match[i] && dfs(i)) {
                    ans++;
                }
            }
        }
        for (int i = 1; i <= n; i++) {
            if (match[i])
                edges.emplace_back(i, match[i] - n);
        }
        return ans;
    }
};

// ==========================================================================================

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    return 0;
}
