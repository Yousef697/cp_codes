#include <bits/stdc++.h>

using namespace std;

// Lowest Common Ancestor (LCA)
struct LCA
{
    int lg = log2(2e5 + 5), cnt = 1;
    vector<int> open, close, level;
    vector<vector<int>> table, adj;

    LCA(int n, int root, vector<vector<int>> &a)
    {
        lg = log2(n);
        open = close = level = vector<int>(n + 1, 0);
        table = vector<vector<int>>(n + 1, vector<int>(lg + 1, root));
        adj = a;

        dfs(root, root);
    }

    void dfs(int node, int parent)
    {
        open[node] = cnt++;

        table[node][0] = parent;
        for (int i = 1; i <= lg; i++)
        {
            int x = table[node][i - 1];
            table[node][i] = table[x][i - 1];
        }

        for (auto i : adj[node])
        {
            if (i != parent)
                level[i] = level[node] + 1, dfs(i, node);
        }

        close[node] = cnt++;
    }

    bool isancestor(int a, int b)
    {
        return open[a] <= open[b] && close[b] <= close[a];
    }

    int query(int a, int b)
    {
        if (isancestor(a, b))
            return a;
        if (isancestor(b, a))
            return b;

        for (int i = lg; i >= 0; i--)
        {
            if (!isancestor(table[a][i], b))
            {
                a = table[a][i];
            }
        }

        return table[a][0];
    }

    int get_distance(int u, int v)
    {
        int ancestor = query(u, v);
        if (ancestor == u || ancestor == v)
            return abs(level[u] - level[v]);
        else
            return abs(level[ancestor] - level[v]) + abs(level[ancestor] - level[u]);
        return 0;
    }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}