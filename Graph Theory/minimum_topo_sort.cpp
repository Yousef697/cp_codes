#include <bits/stdc++.h>
 
using namespace std;
using ll = long long;
 
int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    int n, m;
    cin >> n >> m;
 
    int cnt = 0;
    vector<int> in(n + 1), ans;
    vector<vector<int>> adj(n + 1);
    for (int i = 1; i <= m; i++) {
        int u, v;
        cin >> u >> v;
        in[v]++;
        adj[u].push_back(v);
    }
 
    priority_queue<int, vector<int>, greater<>> pq;
    for (int i = 1; i <= n; i++) {
        if (in[i] == 0) pq.push(i);
    }
 
    while (!pq.empty()) {
        int u = pq.top();
        pq.pop();
 
        ans.push_back(u);
        for (auto v : adj[u]) {
            in[v]--;
            if (in[v] == 0)
                pq.push(v);
        }
        cnt++;
    }
 
    if (cnt != n) {
        cout << "Sandro fails.\n";
    }
    else {
        for (auto i : ans)
            cout << i << " ";
    }
    
    return 0;
}
 