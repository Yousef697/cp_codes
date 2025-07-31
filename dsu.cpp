#include <bits/stdc++.h>

using namespace std;

// Edge
struct edge {
    int u, v, c;

    /// @brief Constructors
    edge() { u = v = c = 0; }
    edge(int a, int b) { u = a, v = b, c = 1; }
    edge(int a, int b, int x) { u = a, v = b, c = x; }

    bool operator<(const edge &e) const { return c < e.c; } // sort edges by weight ascendingly
};

// Union Find, Disjoint Sets Union, DSU
struct UnionFind {
    /// @brief Tree Variables
    int forests, N = 2e5;
    vector<int> parent, size, rank;

    UnionFind() {
        forests = N;
        parent = size = rank = vector<int>(N + 1, 0);
        for (int i = 0; i <= N; i++)
            parent[i] = i, size[i] = rank[i] = 1;
    }

    /// @brief Constructors
    UnionFind(int n) {
        forests = n;
        parent = size = rank = vector<int>(n + 1, 0);
        for (int i = 0; i <= n; i++)
            parent[i] = i, size[i] = rank[i] = 1;
    }

    /// @brief find the root of the tree containing v
    /// @param v
    /// @return the root of the tree containing v
    int find_set(int v) {
        if (v == parent[v])
            return v;
        return parent[v] = find_set(parent[v]);
    }

    /// @brief Union the tree contianing node a with the set containing node b
    /// @param a
    /// @param b
    void union_sets(int a, int b) {
        a = find_set(a);
        b = find_set(b);

        if (a == b)
            return;

        if (size[a] > size[b])
            swap(a, b);

        forests--;
        parent[a] = b;
        size[b] += size[a];
    }

    /// @brief get the number of node of the tree containing node u
    /// @param u
    /// @return size of the tree containing node u
    int size_set(int u) { return size[find_set(u)]; }

    /// @brief check if node u and node v are in the same tree or not
    /// @param u
    /// @param v
    /// @return true if the two nodes are in the same tree, false otherwise
    bool same_set(int u, int v) { return find_set(u) == find_set(v); }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
