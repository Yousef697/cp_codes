#include <bits/stdc++.h>

using namespace std;

struct HLD {
    int n, chain_counter, flat_index;
    vector<int> tree_parent; // parent of node i in the tree
    vector<int> nodes; // number of node in the subtree of node i
    vector<int> leader; // leader of the chain that contain node i
    vector<int> chain_number; // chain number of each node
    vector<int> index_in_chain; // index of each node in its chain
    vector<int> index_in_flat; // index of each node in the flattened tree
    vector<int> flat; // flattened tree
    vector<vector<int>> adj; // adjacency list of each node
    vector<vector<int>> chains; // nodes of each chain

    HLD(const vector<vector<int>>& _adj) {

        adj = _adj;
        flat_index = 0;
        chain_counter = 0;
        n = adj.size() - 1;

        tree_parent = nodes = leader = chain_number = index_in_chain = index_in_flat =
            vector<int> (n + 1);
        chains.resize(n + 1);

        leader[1] = 1;

        dfs_nodes(1, 0);
        dfs_leader(1, 0);
        get_chains(1, 0, 1);

        flat.push_back(0);
        for (int i = 1; i <= n; i++) {
            if (chains[i].empty())
                continue;

            chain_counter++;
            for (auto v : chains[i]) {
                flat.push_back(v);
                index_in_flat[v] = ++flat_index;
            }
        }
    }

    void dfs_nodes(int u, int p) {

        nodes[u] = 1;
        for (auto v : adj[u]) {
            if (v == p)
                continue;

            dfs_nodes(v, u);
            tree_parent[v] = u;
            nodes[u] += nodes[v];
        }
    }

    void dfs_leader(int u, int p) {

        int mx = 0, cnt = 0;
        for (auto v : adj[u])
            if (v != p)
                mx = max(mx, nodes[v]);

        for (auto v : adj[u]) {
            if (v == p)
                continue;

            if (nodes[v] == mx && !cnt)
                leader[v] = leader[u], cnt = 1;
            else
                leader[v] = v;
            dfs_leader(v, u);
        }
    }

    void get_chains(int u, int p, int index) {

        index_in_chain[u] = index;
        chains[leader[u]].push_back(u);

        for (auto v : adj[u]) {
            if (v == p)
                continue;

            if (leader[v] == leader[u])
                get_chains(v, u, index + 1);
            else
                get_chains(v, u, 1);
            /// get_chains(v, u, (leader[v] == leader[u]) * index + 1);
        }
    }

    int query(int u, int v) {
        return 0;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n;
    cin >> n;

    vector<vector<int>> adj(n + 1);

    for (int i = 1; i <= n - 1; i++) {
        int u, v;
        cin >> u >> v;

        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    HLD hld(adj);

    for (int i = 1; i <= n; i++) {
        if (hld.chains[i].empty())
            continue;

        for (auto v : hld.chains[i]) {
            cout << v << ", ";
        }
        cout << "\n";
    }

    for (auto i : hld.flat)
        cout << i << " ";
    cout << "\n";
    for (auto i : hld.index_in_flat)
        cout << i << " ";

    return 0;
}

