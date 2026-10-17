#include <bits/stdc++.h>

using namespace std;
using ll = long long;

const int K = 26;
struct Trie {
    struct Node {
        int cnt = 0;
        vector<int> index;
        array<int, K> next;
    };
    vector<Node> tree;

    Trie() { tree.emplace_back(); }

    int get(char ch) { return ch - 'a'; }
    int advance(int v, char ch) { return tree[v].next[get(ch)]; }

    void insert(const string &s, int idx = 0) {
        int v = 0;
        for (auto& ch : s) {
            tree[v].cnt++;
            if (!advance(v, ch)) {
                tree[v].next[get(ch)] = tree.size();
                tree.emplace_back();
            }
            v = advance(v, ch);
        }
        tree[v].cnt++;
        tree[v].index.emplace_back(idx);
    }
    void erase(const string &s) {
        int v = 0;
        for (auto& ch : s) {
            assert(advance(v, ch));
            tree[v].cnt--;
            v = advance(v, ch);
        }
        tree[v].index.pop_back();
    }
    int search_word(const string &s) {
        int v = 0;
        for (auto& ch : s) {
            v = advance(v, ch);
            if (!v) return 0;
        }
        return tree[v].index.size();
    }
    int search_prefix(const string &s) {
        int v = 0;
        for (auto& ch : s) {
            v = advance(v, ch);
            if (!v) return 0;
        }
        return tree[v].cnt;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
