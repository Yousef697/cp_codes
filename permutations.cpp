#include <bits/stdc++.h>

using namespace std;

// Permutation Cycles
vector<vector<int> > cycles(vector<int> &perm) {
    // 1-indexed
    int n = perm.size() - 1, d;
    vector<int> vis(n + 1);
    vector<vector<int> > ret;
    for (int i = 1; i <= n; i++) {
        if (vis[perm[i]])
            continue;

        vector<int> cycle;
        cycle.push_back(i), vis[i] = 1, d = perm[i];
        while (!vis[d])
            cycle.push_back(d), vis[d] = 1, d = perm[d];
        ret.push_back(cycle);
    }
    return ret;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
