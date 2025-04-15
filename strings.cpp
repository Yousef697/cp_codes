#include <bits/stdc++.h>

using namespace std;

// KMP, String Matching
/*
    prefix: any sub string that starts from index 0
    suffix: any sub string that ends from index size-1
    proper prefix: any prefix that its size is not equal to the main string size

    compute_prefixes function (failure function) compute the proper prefixes of a string such that
    ret[i] = max(proper prefix: prroper prefix = suffix )
*/
vector<int> compute_prefixes(vector<int> pat)
{

    int m = pat.size();

    int i;
    int k;

    vector<int> prefixes(m);

    for (i = 1, k = 0; i < m; i++)
    {

        while (k > 0 && pat[k] != pat[i])
        {
            k = prefixes[k - 1]; // return to the promising position
        }

        if (pat[i] == pat[k])
            prefixes[i] = ++k;
        else
            prefixes[i] = k;
    }

    return prefixes;
}
int KMP(vector<int> s, vector<int> pat)
{

    int n = s.size(), m = pat.size(), ans = 0;

    vector<int> prefixes = compute_prefixes(pat);

    int i; // for s
    int k; // for pat

    for (i = 0, k = 0; i < n; i++)
    {

        while (k > 0 && pat[k] != s[i])
        {
            k = prefixes[k - 1]; // back to the promising position
        }

        if (s[i] == pat[k])
            k++;

        if (k == m)
            ans++, k = prefixes[k - 1]; // make k fail
    }

    return ans;
}

// ==========================================================================================

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}