#include <bits/stdc++.h>

using namespace std;
using ll = long long;

pair<vector<int>, vector<int>> manacher(const string &s) {
    string t = "";
    for (const char &c: s)
        t += '#', t += c;
    t += '#';
    int n = t.size();
    t = '!' + t + '*';

    int l = 0, r = 1;
    vector<int> p(n + 2);
    for (int i = 1; i <= n; i++) {
        p[i] = min(r - i, p[l + (r - i)]);
        p[i] = max(0, p[i]);

        while (t[i - p[i]] == t[i + p[i]])
            p[i]++;

        if (i + p[i] > r)
            l = i - p[i], r = i + p[i];
    }

    for (int &i: p)
        i--;

    auto ret = vector<int>(p.begin() + 1, p.end() - 3);
    vector<int> even, odd;
    for (int i = 0; i + 1 < ret.size(); i += 2) {
        even.push_back(ret[i]);
        odd.push_back(ret[i + 1]);
    }

    // even[i] => if I consider indecies i-1 and i as the center, what is the maximum length of a palindrome substring?
    // odd[i] => if I consider index i as the center, what is the maximum length of a palindrome substring?
    return {even, odd};
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
