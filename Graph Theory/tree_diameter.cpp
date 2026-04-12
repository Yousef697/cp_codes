#include <bits/stdc++.h>

using namespace std;
using ll = long long;

// get tree diameter, assuming the tree is connected
int get_diameter(const vector<vector<int>>& adj) {
    int n = adj.size() - 1;

    auto bfs = [&](int node) -> pair<int, vector<int>> {
        vector<int> dist(n + 1, -1);
        queue<int> q;

        q.push(node);
        dist[node] = 0;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (auto& v : adj[u]) {
                if (dist[v] == -1) {
                    q.push(v);
                    dist[v] = dist[u] + 1;
                }
            }
        }

        int mx = -1, ret = -1;
        for (int i = 1; i <= n; i++) {
            if (dist[i] > mx) {
                ret = i;
                mx = dist[i];
            }
        }
        return {ret, dist};
    };

    auto [_, d] = bfs(1);
    auto [a, da] = bfs(_);
    auto [b, db] = bfs(a);
    // now you have diameter nodes [a, b] with their distance list, do
    return *max_element(da.begin(), da.end());
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
