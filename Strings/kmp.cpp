#include <bits/stdc++.h>

using namespace std;
using ll = long long;

// KMP, String Matching
/*
    prefix: any sub string that starts from index 0
    suffix: any sub string that ends from index size-1
    proper prefix: any prefix that its size is not equal to the main string size

    compute_prefixes function (failure function) compute the proper prefixes of a string such that
    ret[i] = max(proper prefix: prroper prefix = suffix )
*/
vector<int> compute_prefixes(const string& pat) {

    int m = pat.size(), i, k;
    vector<int> prefixes(m);

    for (i = 1, k = 0; i < m; i++) {
        while (k > 0 && pat[k] != pat[i])
            k = prefixes[k - 1]; // return to the promising position

        if (pat[i] == pat[k])
            prefixes[i] = ++k;
        else
            prefixes[i] = k;
    }

    return prefixes;
}
vector<vector<int>> compute_automation(string s) {
    s += '#';
    int n = s.size();
    auto pi = compute_prefixes(s);
    vector<vector<int>> automation(n, vector<int>(26, 0));

    for (int i = 0; i < n; i++) {
        for (int ch = 0; ch < 26; ch++) {
            if (i > 0 && 'a' + ch != s[i])
                automation[i][ch] = automation[pi[i - 1]][ch];
            else
                automation[i][ch] = i + ('a' + ch == s[i]);
        }
    }
    return automation;
}
int KMP(const string& s, const string& pat) {

    int n = s.size(), m = pat.size(), ans = 0;
    vector<int> prefixes = compute_prefixes(pat);

    int i; // for s
    int k; // for pat

    for (i = 0, k = 0; i < n; i++) {
        while (k > 0 && pat[k] != s[i])
            k = prefixes[k - 1]; // back to the promising position

        if (s[i] == pat[k])
            k++;
        if (k == m)
            ans++, k = prefixes[k - 1]; // make k fail
    }

    return ans;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    string s, t;
    cin >> s >> t;

    cout << KMP(s, t);
    
    return 0;
}
