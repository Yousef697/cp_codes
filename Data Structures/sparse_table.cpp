#include <bits/stdc++.h>

using namespace std;

// Sparse Table
struct SparseTable {
    int n;
    vector<int> logs;
    vector<vector<int> > table;

    int merge(int a, int b) { return __gcd(a, b); }

    SparseTable() {}

    SparseTable(vector<int> &vec) {
        n = vec.size() - 1;
        logs.assign(n + 1, 0);

        logs[1] = 0;
        for (int i = 2; i <= n; i++)
            logs[i] = logs[i / 2] + 1;
        table = vector<vector<int> >(n + 1, vector<int>(logs[n] + 2, 0));

        for (int i = 1; i <= n; i++)
            table[i][0] = vec[i];
        for (int j = 1; j <= logs[n]; j++)
            for (int i = 1; i <= n - (1 << j) + 1; i++)
                table[i][j] = merge(table[i][j - 1], table[i + (1 << (j - 1))][j - 1]);
    }

    // Min, Max, GCD (we can take the element more than once)
    int fast_query(int l, int r) {
        int len = r - l + 1;
        return merge(table[l][logs[len]], table[r - (1 << logs[len]) + 1][logs[len]]);
    }

    // Disjoint Ranges: Sums, XORs (we cannot take the element more than once)
    int slow_query(int l, int r) {

        int len = r - l + 1, k = 0, ans = 0;
        while (len) {
            if (len & 1)
                ans = merge(ans, table[l][k]), l += (1 << k);
            k++, len >>= 1;
        }
        return ans;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
