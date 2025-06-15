#include <bits/stdc++.h>

using namespace std;

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n;
    cin >> n;

    vector<int> tree_parent(n + 1); // parent of node i in the tree
    vector<vector<int>> adj(n + 1); // adjacency list of each node
    for (int i = 1; i <= n - 1; i++) {
        int u, v;
        cin >> u >> v;

        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector<int> dp_nodes(n + 1); // number of node in the subtree of node i
    function<void(int, int)> dfs_nodes = [&](int u, int p) {

        dp_nodes[u] = 1;
        for (auto v : adj[u]) {
            if (v == p)
                continue;

            dfs_nodes(v, u);
            tree_parent[v] = u;
            dp_nodes[u] += dp_nodes[v];
        }
    };
    dfs_nodes(1, 0);

    vector<int> leader(n + 1); leader[1] = 1; // leader on the chain that contain node i
    function<void(int, int)> dfs_leader = [&](int u, int p) {

        int mx = 0, cnt = 0;
        for (auto v : adj[u])
            if (v != p)
                mx = max(mx, dp_nodes[v]);

        for (auto v : adj[u]) {
            if (v == p)
                continue;

            if (dp_nodes[v] == mx && !cnt)
                leader[v] = leader[u], cnt = 1;
            else
                leader[v] = v;
            dfs_leader(v, u);
        }
    };
    dfs_leader(1, 0);

    int chain_counter = 1;
    vector<vector<int>> chains(n + 1);
    vector<int> chain_number(n + 1), rank_in_chain(n + 1);
    function<void(int, int, int)> get_chains = [&](int u, int p, int rank) {

        rank_in_chain[u] = rank;
        chains[leader[u]].push_back(u);

        for (auto v : adj[u]) {
            if (v == p)
                continue;

            if (leader[v] == leader[u])
                get_chains(v, u, rank + 1);
            else
                get_chains(v, u, 1), chain_counter++;
            /// get_chains(v, u, (leader[v] == leader[u]) * rank + 1);
        }
    };
    get_chains(1, 0, 1);

    // printing the chains
    int j = 0;
    for (int i = 1; i <= n; i++) {
        if (chains[i].empty())
            continue;

        cout << "Chain " << ++j << ": ";
        for (auto v : chains[i])
            cout << v << ", ";
        cout << "\n";
    }

    return 0;
}
