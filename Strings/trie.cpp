#include <bits/stdc++.h>

using namespace std;
using ll = long long;

const int K = 2;
struct Trie {
    vector<int> leaf;
    vector<int> cnts;
    vector<vector<int> > tree;

    Trie() {
        leaf.push_back(0);
        cnts.push_back(0);
        tree.push_back(vector<int>(K, 0));
    }

    void insert(const string &s, int v = 0, int i = 0) {
        if (i == s.size()) {
            cnts[v]++;
            leaf[v] = 1;
            return;
        }

        if (!tree[v][s[i] - '0']) {

            tree[v][s[i] - '0'] = tree.size();
            leaf.push_back(0);
            cnts.push_back(0);
            tree.push_back(vector<int>(2, 0));
        }

        int next = tree[v][s[i] - '0'];
        insert(s, next, i + 1);
        cnts[v]++;
    }
    void erase(const string &s, int v = 0, int i = 0) {
        if (i == s.size()) {
            cnts[v]--;
            if (cnts[v] == 0)
                leaf[v] = 0;
            return;
        }

        assert(tree[v][s[i] - '0']);
        int next = tree[v][s[i] - '0'];
        erase(s, next, i + 1);

        cnts[v]--;
        if (cnts[next] == 0)
            tree[v][s[i] - '0'] = 0;
    }
    int serach_word(const string &s, int v = 0, int i = 0) {
        if (i == s.size())
            return leaf[v];

        if (!tree[v][s[i] - '0'])
            return -1;
        return serach_word(s, tree[v][s[i] - '0'], i + 1);
    }
    int serach_prefix(const string &s, int v = 0, int i = 0) {
        if (i == s.size())
            return 1;

        if (!tree[v][s[i] - '0'])
            return -1;
        return serach_prefix(s, tree[v][s[i] - '0'], i + 1);
    }
};


int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
