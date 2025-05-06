#include <bits/stdc++.h>

using namespace std;

struct fenwick
{
    int n;
    vector<int> bit;

    fenwick(int _n)
    {
        n = _n;
        bit = vector<int> (n + 1, 0);
    }

    fenwick(const vector<int>& a)
    {
        n = a.size() - 1;
        bit = vector<int> (n + 1, 0);

        for (int i = 1; i <= n; i++)
        {
            bit[i] = bit[i] + a[i];
            int next = i + (i & (-i));
            if (next <= n)
                bit[next] = bit[next] + bit[i];
        }
    }

    void update(int idx, int val)
    {
        while (idx <= n)
        {
            bit[idx] = bit[idx] + val;
            idx += ((+idx) & (-idx));
        }
    }

    int query(int idx)
    {
        int ans = 0;
        while (idx)
        {
            ans = ans + bit[idx];
            idx -= ((+idx) & (-idx));
        }
        return ans;
    }
};

struct fenwick2d
{
    int n, m;
    vector<vector<int>> bit;

    fenwick2d(int _n, int _m)
    {
        n = _n, m = _m;
        bit = vector<vector<int>> (n + 1, vector<int> (m + 1, 0));
    }

    void update(int idxX, int idxY, int val)
    {
        while (idxX <= n)
        {
            int tmp = idxY;
            while (tmp <= n)
            {
                bit[idxX][tmp] += val;
                tmp += ((+tmp) & (-tmp));
            }
            idxX += ((+idxX) & (-idxX));
        }
    }

    int query(int idxX, int idxY)
    {
        int ans = 0;
        while (idxX)
        {
            int tmp = idxY;
            while (tmp)
            {
                ans += bit[idxX][tmp];
                tmp -= ((+tmp) & (-tmp));
            }
            idxX -= ((+idxX) & (-idxX));
        }
        return ans;
    }
};

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    
    
    return 0;
}
