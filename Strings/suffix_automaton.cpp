#include <bits/stdc++.h>

using namespace std;
using ll = long long;

struct SuffixAutomaton {

    const int ALPHA = 256;
    int lst, sz;
    vector<int> len, fail;
    vector<map<char, int>> child;

    SuffixAutomaton(const string& s) {
        build(s);
    }

    int add_node() {
        map<char, int> temp;
        child.push_back(temp);
        len.emplace_back(0), fail.emplace_back(0);
        return sz++;
    }
    void init() {
        len.resize(2);
        fail.resize(2);
        child.resize(2);

        lst = 1, sz = 2;
        len[1] = 1;
        fail[1] = len[0] = 0;

        for (int i = 0; i < ALPHA; i++) {
            child[0][i] = 1;
        }
    }
    void build(const string& s) {
        init();
        for (auto i : s) {
            add_char(i);
            // print();
        }
    }

    void add_char(unsigned char c) {
        int u = add_node();
        len[u] = len[lst] + 1;

        int v, w;
        for (v = lst;; v = fail[v]) {
            w = child[v].emplace(c, u).first->second;
            if (w != u) break;
        }
        if (len[v] + 1 == len[w]) {
            fail[u] = w;
        }
        else {
            int x = add_node();
            len[x] = len[v] + 1;
            child[x] = child[w];
            int y;
            for (y = v;; y = fail[y]) {
                auto& z = child[y][c];
                if (z == w) z = x;
                else break;
            }
            fail[x] = fail[w];
            fail[w] = fail[u] = x;
        }
        lst = u;
    }

    void print() {
        for (int u = 1; u < sz; u++) {
            for (auto& [c, v] : child[u]) {
                cout << u << "," << len[u] << "," << fail[u] << " ";
                cout << v << "," << len[v] << "," << fail[v] << " ";
                cout << c << "\n";
            }
        }
        cout << "------------------------------------------------------\n";
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
