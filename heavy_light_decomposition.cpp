#include <bits/stdc++.h>

using namespace std;

template <typename T>
struct SegmentTree
{
    int size;
    vector<T> v, seg;

    SegmentTree() {}

    SegmentTree(vector<T> vec)
    {
        v = vec;
        size = 1;
        int n = vec.size() - 1;
        while (size < n)
            size *= 2;
        seg.assign(2 * size + 2, 0);
        build();
    }

    T merge(T a, T b)
    {
        return min(a, b);
    }

    void build(int n, int l, int r)
    {
        if (l == r)
        {
            if (l < v.size())
                seg[n] = v[l];
            return;
        }

        build(2 * n, l, (l + r) / 2);
        build(2 * n + 1, (l + r) / 2 + 1, r);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    void update(int n, int l, int r, int ind, int val)
    {
        if (ind < l || ind > r)
            return;

        if (l == r)
        {
            if (l < v.size())
                seg[n] = val, v[l] = val;
            return;
        }

        update(2 * n, l, (l + r) / 2, ind, val);
        update(2 * n + 1, (l + r) / 2 + 1, r, ind, val);

        seg[n] = merge(seg[2 * n], seg[2 * n + 1]);
    }

    T query(int n, int l, int r, int lq, int rq)
    {
        if (r < lq || l > rq)
            return INT_MAX;

        if (lq <= l && r <= rq)
            return seg[n];

        T a = query(2 * n, l, (l + r) / 2, lq, rq);
        T b = query(2 * n + 1, (l + r) / 2 + 1, r, lq, rq);

        return merge(a, b);
    }

    void build() { build(1, 1, size); }

    void update(int ind, int val) { update(1, 1, size, ind, val); }

    T query(int l, int r) { return query(1, 1, size, l, r); }
};

struct HLD {
    int n, chain_counter, flat_index;
    vector<int> tree_parent; // parent of node i in the tree
    vector<int> depth; // depth of each node
    vector<int> nodes; // number of node in the subtree of node i
    vector<int> leader; // leader of the chain that contain node i
    vector<int> chain_number; // chain number of each node
    vector<int> index_in_chain; // index of each node in its chain
    vector<int> index_in_flat; // index of each node in the flattened tree
    vector<int> flat; // flattened tree
    vector<int> values; // value of node at index flat[i]
    vector<vector<int>> adj; // adjacency list of each node
    vector<vector<int>> chains; // nodes of each chain
    SegmentTree<int> tree;

    HLD(const vector<vector<int>>& _adj, const vector<int>& vals) {

        adj = _adj;
        flat_index = 0;
        chain_counter = 0;
        n = adj.size() - 1;

        tree_parent = depth = nodes = leader = chain_number = index_in_chain = index_in_flat =
            vector<int> (n + 1);
        chains.resize(n + 1);

        depth[1] = 0;
        leader[1] = 1;

        dfs_nodes(1, 0);
        dfs_leader(1, 0);
        get_chains(1, 0, 1);

        flat.push_back(0);
        values.push_back(0);

        for (int i = 1; i <= n; i++) {
            if (chains[i].empty())
                continue;

            chain_counter++;
            for (auto v : chains[i]) {
                flat.push_back(v);
                values.push_back(vals[v]);
                index_in_flat[v] = ++flat_index;
            }
        }
        tree = SegmentTree<int>(values);
    }

    void dfs_nodes(int u, int p) {

        nodes[u] = 1;
        for (auto v : adj[u]) {
            if (v == p)
                continue;

            dfs_nodes(v, u);
            tree_parent[v] = u;
            nodes[u] += nodes[v];
            depth[v] = depth[u] + 1;
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

    void update(int u, int val) {
        int ind = index_in_flat[u];
        tree.update(ind, val);
    }

    int query(int u, int v) {

        int res = INT_MAX, l, r;

        while (leader[u] != leader[v]) {
            if (depth[u] > depth[v])
                swap(u, v);

            l = index_in_flat[leader[v]];
            r = index_in_flat[v];
            res = min(res, tree.query(l, r));
            v = tree_parent[leader[v]];
        }

        if (depth[u] > depth[v])
            swap(u, v);
        l = index_in_flat[u];
        r = index_in_flat[v];
        if (l > r)
            swap(l, r);
        res = min(res, tree.query(l, r));

        return res;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n;
    cin >> n;

    vector<int> values(n + 1);
    vector<vector<int>> adj(n + 1);

    for (int i = 1; i <= n; i++) {
        cin >> values[i];
    }

    for (int i = 1; i <= n - 1; i++) {
        int u, v;
        cin >> u >> v;

        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    HLD hld(adj, values);

    /*
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
    */

    int q;
    cin >> q;

    while (q--) {
        int type, u, v;
        cin >> type >> u >> v;

        if (type == 1) {
            hld.update(u, v);
        }
        else {
            cout << hld.query(u, v) << "\n";
        }
    }

    return 0;
}

/*
20
4 5 3 6 8 3 5 6 3 5 6 8 1 9 8 6 2 3 5 3
1 2
1 3
1 4
2 5
3 8
3 9
4 10
5 6
5 7
10 11
10 12
10 13
11 14
12 15
14 16
14 17
14 18
7 19
7 20
11
2 7 9
2 8 10
2 6 8
2 7 9
2 10 5
2 20 2
2 13 15
1 3 1
2 8 9
2 8 6
2 9 1

01 02 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20
4  5  3  6  8  3  5  6  3  5  6  8  1  9  8  6  2  3  5  3

*/