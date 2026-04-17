#include <bits/stdc++.h>
#define int long long

using namespace std;
using ll = long long;
const int N = 2e5 + 5;

int n;
int subtree[N];
int down[N], up[N], answer[N];
vector<int> adj[N];

int merge(int a, int b) {
    return (a + b);
}
void dfs_down(int u, int p) {
    down[u] = 0;
    subtree[u] = 1;

    for (auto& v : adj[u]) {
        if (v == p) continue;
        dfs_down(v, u);
        subtree[u] += subtree[v];
        down[u] = merge(down[u], down[v] + subtree[v]); // this logic may differ
    }
}
void dfs_up(int u, int p) {
    vector<int> prefix, suffix;

    // this logic may differ
    for (auto& v : adj[u]) {
        int val = 0;
        if (v != p) val = down[v] + 2 * subtree[v];
        else val = up[u] + n - subtree[u] + 1;
        prefix.push_back(val);
        suffix.push_back(val);
    }

    int sz = prefix.size();
    for (int i = 1; i < sz; i++) prefix[i] = merge(prefix[i], prefix[i - 1]);
    for (int i = sz - 2; i >= 0; i--) suffix[i] = merge(suffix[i], suffix[i + 1]);

    // this logic may differ
    for (int i = 0; i < sz; i++) {
        int v = adj[u][i];
        if (v == p) continue;
        int val1 = (i ? prefix[i - 1] : 0);
        int val2 = (i != sz - 1 ? suffix[i + 1] : 0);
        up[v] = merge(val1, val2);
    }

    for (auto& v : adj[u]) {
        if (v == p) continue;
        dfs_up(v, u);
    }
}
void reroot() {
    dfs_down(1, 0);
    dfs_up(1, 0);
    for (int i = 1; i <= n; i++) {
        answer[i] = merge(down[i], up[i]); // f(down[i], up[i])
        if (i != 1) answer[i]++;
    }
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    cin >> n;
    for (int i = 1; i < n; i++) {
        int u, v;
        cin >> u >> v;
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    reroot();

    for (int i = 1; i <= n; i++) {
        // cout << i << ": " << answer[i] << "\n";
        cout << answer[i] << " ";
    }
    
    return 0;
}
