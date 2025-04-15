#include <bits/stdc++.h>

using namespace std;

// Sqrt Decomposition
// For Sum
struct SqrtDecomposition
{
    int n, sq = 448;
    vector<long long> arr, ans;
    vector<vector<long long>> blocks;

    SqrtDecomposition(vector<long long> &v)
    {
        arr = v;
        n = v.size();
        ans.assign(sq, 0);
        blocks.resize(sq);

        for (int i = 0; i < n; i++)
        {
            ans[i / sq] += arr[i];
            blocks[i / sq].emplace_back(arr[i]);
        }
    }

    void update(int idx, long long val)
    {
        int i = idx / sq, j = idx % sq;

        ans[i] -= arr[idx];
        arr[idx] = val;
        ans[i] += arr[idx];
        blocks[i][j] = val;
    }

    long long query(int l, int r)
    {
        long long sum = 0;

        while (l <= r)
        {
            if (l % sq == 0 && l + sq - 1 <= r)
                sum += ans[l / sq], l += sq;
            else
                sum += arr[l], l++;
        }

        return sum;
    }
};

// For Min
struct SqrtDecomposition
{
    int n, sq = 448;
    vector<int> arr, ans;
    vector<vector<int>> blocks;

    SqrtDecomposition(vector<int> &v)
    {
        arr = v;
        n = v.size();
        ans.assign(sq, 1e9);
        blocks.resize(sq);

        for (int i = 0; i < n; i++)
        {
            ans[i / sq] = min(ans[i / sq], arr[i]);
            blocks[i / sq].emplace_back(arr[i]);
        }
    }

    void update(int idx, int val)
    {
        int i = idx / sq, j = idx % sq, x = 1e9;

        arr[idx] = val;
        blocks[i][j] = val;

        for (auto k : blocks[i])
            x = min(x, k);

        ans[idx / sq] = x;
    }

    int query(int l, int r)
    {
        int ret = 1e9;

        while (l <= r)
        {
            if (l % sq == 0 && l + sq - 1 <= r)
                ret = min(ret, ans[l / sq]), l += sq;
            else
                ret = min(ret, arr[l]), l++;
        }

        return ret;
    }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}