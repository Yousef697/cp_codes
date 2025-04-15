#include <bits/stdc++.h>

using namespace std;

// 2-sat, 2satisfiability
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
struct TwoSat
{
    vector<int> s;
    int n, cnt = 1;
    vector<vector<int>> adj, scc;
    vector<int> id, low, visited;

    TwoSat(int v, vector<vector<int>> &a)
    {
        // node i: Xi, node i + n: !Xi
        n = v, adj = a;
        id = low = visited = vector<int>(n + 1, 0);

        for (int i = 1; i <= n; i++)
        {
            if (id[i] == 0)
                dfs(i);
        }
    }

    void dfs(int node)
    {
        s.push_back(node);
        visited[node] = 1;
        id[node] = low[node] = cnt++;

        for (int neighbor : adj[node])
        {
            if (id[neighbor] == 0)
                dfs(neighbor);
            if (visited[neighbor])
                low[node] = min(low[node], low[neighbor]);
        }

        if (id[node] == low[node])
        {
            vector<int> component;
            while (s.size())
            {
                int x = s.back();
                component.push_back(x);
                visited[x] = 0;
                s.pop_back();

                if (x == node)
                    break;
            }
            scc.push_back(component);
        }
    }

    vector<int> solve()
    {
        UnionFind uf(n);
        for (vector<int> &component : scc)
        {
            for (int u : component)
                uf.union_sets(component[0], u);
        }

        int ok = 1;
        for (int i = 1; i <= n / 2; i++)
            if (uf.same_set(i, i + n / 2))
                ok = 0;
        if (!ok)
            return vector<int>(n / 2 + 1, -1);

        vector<vector<int>> new_adj(n + 1), components(n + 1);

        for (int u = 1; u <= n; u++)
        {
            for (int v : adj[u])
                if (!uf.same_set(u, v))
                    new_adj[uf.find_set(u)].push_back(uf.find_set(v));
        }

        for (vector<int> &component : scc)
        {
            for (int u : component)
                components[uf.find_set(u)].push_back(u);
        }

        vector<int> vis(n + 1), topo;

        function<void(int)> toposort = [&](int u)
        {
            vis[u] = 1;
            for (int v : new_adj[u])
            {
                if (!vis[v])
                    toposort(v);
            }
            topo.push_back(u);
        };

        for (int i = 1; i <= n; i++)
            if (!vis[i])
                toposort(i);

        vector<int> ans(n / 2 + 1, -1);
        for (int u : topo)
        {
            for (int v : components[u])
            {
                if (v <= n / 2 && ~ans[v])
                    continue;
                if (v > n / 2 && ~ans[v - n / 2])
                    continue;

                if (v <= n / 2)
                    ans[v] = 1;
                else
                    ans[v - n / 2] = 0;
            }
        }

        return ans;
    }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, m; // edges, nodes
    cin >> n >> m;

    vector<vector<int>> adj(2 * m + 1);
    for (int i = 1; i <= n; i++)
    {
        char a, b;
        int u, v, nu, nv; // u, v, !u, !v
        cin >> a >> u >> b >> v;

        nu = u + m, nv = v + m;
        if (a == '-')
            swap(u, nu);
        if (b == '-')
            swap(v, nv);

        adj[nv].push_back(u); // edge from !v to u
        adj[nu].push_back(v); // edges from !u to v
    }

    TwoSat ret(2 * m, adj);
    vector<int> ans = ret.solve();

    if (ans[1] == -1)
        cout << "IMPOSSIBLE";
    else
        for (int i = 1; i <= m; i++)
            cout << (ans[i] ? '+' : '-') << " ";

    return 0;
}