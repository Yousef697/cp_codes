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

// Dinic
struct Dinic {
    const int inf = 1e18 + 18;
    struct Edge {
        int u, v, cap, flow = 0;
        Edge(int u, int v, int cap): u(u), v(v), cap(cap) {}
    };

    int n, m = 0;
    vector<Edge> edges, init;
    vector<vector<int>> adj;
    vector<int> level, ptr;
    queue<int> q;

    Dinic(int n): n(n) {
        adj = vector<vector<int>>(n + 1);
        ptr = level = vector<int>(n + 1);
    }

    void add_edge(int u, int v, int cap) {
        init.emplace_back(u, v, cap);
        init.emplace_back(v, u, 0);
        adj[u].emplace_back(m); m++;
        adj[v].emplace_back(m); m++;
    }
    bool bfs(int s, int t) {
        level = vector<int>(n + 1, -1);
        level[s] = 0;
        q.emplace(s);
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (auto& id : adj[u]) {
                if (edges[id].cap == edges[id].flow) continue;
                if (level[edges[id].v] != -1) continue;
                level[edges[id].v] = level[u] + 1;
                q.emplace(edges[id].v);
            }
        }
        return level[t] != -1;
    }
    int dfs(int u, int pushed, int t) {
        if (pushed == 0) return 0;
        if (u == t) return pushed;
        for (int& cid = ptr[u]; cid < adj[u].size(); cid++) {
            int id = adj[u][cid];
            int v = edges[id].v;
            if (level[u] + 1 != level[v]) continue;
            int new_pushed = min(pushed, edges[id].cap - edges[id].flow);
            int add = dfs(v, new_pushed, t);
            if (add == 0) continue;
            edges[id].flow += add;
            edges[id ^ 1].flow -= add;
            return add;
        }
        return 0;
    }
    int max_flow(int s, int t) {
        int ans = 0;
        edges = init;

        while (true) {
            if (!bfs(s, t)) break;
            ptr = vector<int>(n + 1, 0);
            while (true) {
                int pushed = dfs(s, inf, t);
                if (!pushed) break;
                ans += pushed;
            }
        }
        return ans;
    }
};

// ==========================================================================================

// Min-Cost Max-Flow
struct MCMF {
    struct edge {
        int to, cap, cost, rev;
    };
    const int inf = 1e18;

    int n;
    vector<vector<edge>> adj;
    vector<int> dist, parent_node, parent_edge;

    MCMF(int n): n(n) {
        adj = vector<vector<edge>>(n + 1);
        dist = parent_node = parent_edge = vector<int>(n + 1);
    }

    void add_edge(int u, int v, int cap, int cost) {
        adj[u].push_back({v, cap, cost, (int)adj[v].size()});
        adj[v].push_back({u, 0, -cost, (int)adj[u].size() - 1});
    }

    pair<int, int> mcmf(int s, int t, int k = -1) {
        if (k == -1) k = inf;
        int flow = 0, cost = 0;
        while (flow < k) {
            fill(dist.begin(), dist.end(), inf);
            priority_queue<pair<int, int>, vector<pair<int, int>>, greater<>> pq;

            dist[s] = 0;
            pq.emplace(0, s);

            while (!pq.empty()) {
                auto [d, u] = pq.top();
                pq.pop();
                if (d > dist[u]) continue;

                for (int i = 0; i < adj[u].size(); i++) {
                    auto& e = adj[u][i];
                    if (e.cap > 0 && dist[e.to] > d + e.cost) {
                        dist[e.to] = d + e.cost;
                        parent_node[e.to] = u;
                        parent_edge[e.to] = i;
                        pq.emplace(dist[e.to], e.to);
                    }
                }
            }
            if (dist[t] == inf) break;
            int add = k - flow, cur = t;
            while (cur != s) {
                auto& e = adj[parent_node[cur]][parent_edge[cur]];
                add = min(add, e.cap);
                cur = parent_node[cur];
            }
            cur = t;
            while (cur != s) {
                auto& e = adj[parent_node[cur]][parent_edge[cur]];
                e.cap -= add;
                adj[cur][e.rev].cap += add;
                cost += add * e.cost;
                cur = parent_node[cur];
            }
            flow += add;
        }

        return {flow, cost};
    }
};

// ==========================================================================================

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    return 0;
}
