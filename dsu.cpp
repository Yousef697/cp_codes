#include <bits/stdc++.h>

using namespace std;

// Edge
struct edge
{
    int u, v, c;

    edge() {}
    edge(int a, int b) { u = a, v = b, c = 1; }
    edge(int a, int b, int x) { u = a, v = b, c = x; }

    bool operator<(const edge &e) { return c < e.c; }
};

// Union Find, Disjoint Sets Union, DSU
struct UnionFind
{
    int forests;
    vector<int> parent, size, rank;

    UnionFind(int n)
    {
        forests = n;
        parent = size = rank = vector<int>(n + 1, 0);
        for (int i = 0; i <= n; i++)
            parent[i] = i, size[i] = rank[i] = 1;
    }

    int find_set(int v)
    {
        if (v == parent[v])
            return v;
        return parent[v] = find_set(parent[v]);
    }

    void union_sets(int a, int b)
    {
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

    int size_set(int u) { return size[find_set(u)]; }

    bool same_set(int u, int v) { return find_set(u) == find_set(v); }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}