#include <bits/stdc++.h>

using namespace std;

// Articulation Points and Bridges
struct ArticulationPointsAndBridges
{
    vector<vector<int>> adj;
    vector<pair<int, int>> bridges;
    int n, m, cnt = 1, root, root_children;
    vector<int> id, low, parent, articulation;

    ArticulationPointsAndBridges(vector<vector<int>> &a)
    {
        adj = a, n = a.size() - 1, cnt = 1;
        id = low = parent = articulation = vector<int>(n + 1, 0);

        for (int i = 1; i <= n; i++)
        {
            if (id[i] == 0)
            {
                root = i, root_children = 0;
                dfs(i);

                articulation[i] = (root_children > 1);
            }
        }
    }

    void dfs(int node)
    {
        id[node] = low[node] = cnt++;

        for (int neighbor : adj[node])
        {
            if (neighbor == parent[node])
                continue;
            if (id[neighbor] == 0)
            {
                parent[neighbor] = node;
                dfs(neighbor);

                if (node == root)
                    root_children++;

                low[node] = min(low[node], low[neighbor]);
                if (low[neighbor] >= id[node])
                    articulation[node] = 1;
                if (low[neighbor] > id[node])
                    bridges.push_back({node, neighbor});
            }
            else
                low[node] = min(low[node], id[neighbor]);
        }
    }

    vector<pair<int, int>> GetBridges() { return bridges; }

    vector<int> GetArticulationPoints()
    {
        vector<int> ret;
        for (int i = 1; i <= n; i++)
            if (articulation[i])
                ret.push_back(i);
        return ret;
    }
};

// ==========================================================================================

// Tarjan, SCC
struct Tarjan
{
    vector<int> s;
    int n, m, cnt = 1;
    vector<vector<int>> adj, scc;
    vector<int> id, low, visited;

    Tarjan(vector<vector<int>> &a)
    {
        n = a.size() - 1, adj = a;
        id = low = visited = vector<int>(n + 1, 0);

        for (int i = 1; i <= n; i++)
        {
            if (id[i] == 0)
            {
                dfs(i);
            }
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

    vector<vector<int>> GetSCC() { return scc; }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    int n, m;
    cin >> n >> m;

    vector<vector<int>> adj(n + 1);
    for (int i = 1; i <= m; i++)
    {
        int u, v;
        cin >> u >> v;

        adj[u].push_back(v);
    }

    Tarjan ret(adj);
    ArticulationPointsAndBridges ret2(adj);

    vector<vector<int>> scc = ret.GetSCC();
    vector<int> points = ret2.GetArticulationPoints();
    vector<pair<int, int>> bridges = ret2.GetBridges();

    // Do what you want on bridges, points, scc

    return 0;
}