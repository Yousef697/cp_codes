#include <bits/stdc++.h>

using namespace std;

// Sqrt Decomposition (Sum)
struct SqrtDecomposition {
    int n, sq = 448;
    vector<int> arr, ans;
    vector<vector<int> > blocks;

    SqrtDecomposition(vector<int> &vec) {
        arr = vec;
        n = vec.size();
        ans = vector<int>(sq, 0);
        blocks = vector<vector<int>>(sq);

        for (int i = 0; i < n; i++) {
            ans[i / sq] += vec[i];
            blocks[i / sq].emplace_back(vec[i]);
        }
    }

    void update(int idx, int val) {
        int i = idx / sq, j = idx % sq;

        ans[i] -= arr[idx];
        arr[idx] = val;
        ans[i] += arr[idx];
        blocks[i][j] = val;
    }
    int query(int l, int r) {

        int sum = 0;
        while (l <= r) {
            if (l % sq == 0 && l + sq - 1 <= r)
                sum += ans[l / sq], l += sq;
            else
                sum += arr[l], l++;
        }

        return sum;
    }
};

// Sqrt Decomposition (Min, Max, Gcd)
struct Sqrt_Decomposition {
    int n, sq = 448;
    vector<int> arr, ans;
    vector<vector<int> > blocks;

    Sqrt_Decomposition(vector<int> &vec) {
        arr = vec;
        n = vec.size();
        ans = vector<int>(sq, 1e9);
        blocks = vector<vector<int>>(sq);

        for (int i = 0; i < n; i++) {
            ans[i / sq] = min(ans[i / sq], arr[i]);
            blocks[i / sq].emplace_back(arr[i]);
        }
    }

    int merge(int a, int b) { return min(a, b); }

    void update(int idx, int val) {

        int i = idx / sq, j = idx % sq, x = 1e9;
        arr[idx] = val;
        blocks[i][j] = val;

        for (auto k: blocks[i])
            x = merge(x, k);

        ans[idx / sq] = x;
    }
    int query(int l, int r) {

        int ret = 1e9;
        while (l <= r) {
            if (l % sq == 0 && l + sq - 1 <= r)
                ret = merge(ret, ans[l / sq]), l += sq;
            else
                ret = merge(ret, arr[l]), l++;
        }

        return ret;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}
