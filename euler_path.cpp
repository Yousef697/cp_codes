#include <bits/stdc++.h>

using namespace std;

// Eulerian Circuit, Euler Tour, Euler Path
struct UndirectedEulerCircuit {
    int n = -1, m = -1, ok = 1;
    int u1 = -1, u2 = -1; // two nodes with odd degree

    vector<pair<int, int> > edges;
    vector<vector<pair<int, int> > > adj;
    vector<int> degree, ans, vis, taken;

    UndirectedEulerCircuit(vector<pair<int, int> > &_edges) {
        edges = _edges;
        m = edges.size() - 1;

        /// get number of nodes
        for (auto [u, v]: edges)
            n = max({n, u, v});

        /// initialize the graph
        adj.resize(n + 1);
        vis = vector<int>(n + 1, 0);
        taken = vector<int>(m + 1, 0);
        degree = vector<int>(n + 1, 0);

        for (int i = 1; i <= m; i++) {
            auto [u, v] = edges[i];
            adj[u].push_back({v, i});
            adj[v].push_back({u, i});
            degree[u]++, degree[v]++;
        }

        for (int i = 1; i <= n; i++) {
            if (degree[i] % 2 == 1) {
                if (u1 == -1)
                    u1 = i;
                else if (u2 == -1)
                    u2 = i;
                else
                    ok = 0;
            }
        }

        if (u1 != -1 && u2 == -1)
            ok = 0;
    }

    void dfs(int u) {
        while (!adj[u].empty()) {
            auto [v, ind] = adj[u].back();
            adj[u].pop_back();

            if (taken[ind])
                continue;

            taken[ind] = 1, dfs(v);
        }
        ans.push_back(u);
    }

    vector<int> solve() {
        if (!ok)
            return vector<int>{-1};

        /// if there is only two node with odd degree, add a virtual edge between these two nodes
        if (u1 != -1)
            m++, adj[u1].push_back({u2, m}), adj[u2].push_back({u1, m});

        taken = vector<int>(m + 1, 0);
        dfs(1);
        for (int i = 1; i <= m; i++)
            if (!taken[i])
                return vector<int>{-1};
        return ans;
    }
};

struct DirectedEulerCircuit {
    int n = -1, m = -1, ok = 1;

    vector<pair<int, int> > edges;
    vector<vector<pair<int, int> > > adj;
    vector<int> in, out, ans, vis, taken;

    DirectedEulerCircuit(vector<pair<int, int> > &_edges) {
        edges = _edges;
        m = edges.size() - 1;

        for (auto [u, v]: edges)
            n = max({n, u, v});

        adj.resize(n + 1);
        in = vector<int>(n + 1, 0);
        out = vector<int>(n + 1, 0);
        vis = vector<int>(n + 1, 0);
        taken = vector<int>(m + 1, 0);

        for (int i = 1; i <= m; i++) {
            auto [u, v] = edges[i];
            adj[u].push_back({v, i});
            out[u]++, in[v]++;
        }

        for (int i = 2; i < n; i++) {
            if (in[i] != out[i])
                ok = 0;
        }
        // circuit from 1 to n (can be changed)
        if (out[1] != in[1] + 1 || in[n] != out[n] + 1)
            ok = 0;
    }

    void dfs(int u) {
        while (!adj[u].empty()) {
            auto [v, ind] = adj[u].back();
            adj[u].pop_back();

            if (taken[ind])
                continue;

            taken[ind] = 1;
            dfs(v);
        }
        ans.push_back(u);
    }

    vector<int> solve() {
        if (!ok)
            return vector<int>{-1};

        taken = vector<int>(m + 1, 0);
        dfs(1);
        for (int i = 1; i <= m; i++)
            if (!taken[i])
                return vector<int>{-1};
        reverse(ans.begin(), ans.end());
        return ans;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
