#include <bits/stdc++.h>

using namespace std;
using ll = long long;

// My Code (Have issues)
const int K = 128, BASE = 'A';
struct vertex {
    char ch;
    bool leaf = false;
    int next[K], parent, fail, exit, id;
    vector<int> ids;

    vertex(int _p = -1, char _ch = '$') {
        ch = _ch;
        parent = _p;
        fail = exit = 0;
        memset(next, -1, sizeof next);
    }
};
struct aho_corasick {
    int time = -1;
    vector<vertex> tree;
    vector<int> in, out;
    aho_corasick() { tree.emplace_back(); }

    void insert(const string &s, int id = 0) {
        int u = 0;
        for (auto &c: s) {
            int v = c - BASE;
            if (tree[u].next[v] == -1) {
                tree[u].next[v] = tree.size();
                tree.emplace_back(u, c);
            }
            u = tree[u].next[v];
        }
        // tree[u].id = id;
        tree[u].ids.emplace_back(id);
        tree[u].leaf = true;
    }
    void build_aho() {
        queue<int> q;
        q.push(0);

        auto compute_failure = [&](int u) -> int {
            if (u == 0 || tree[u].parent == 0)
                return 0;

            int x = tree[u].ch - BASE;
            while (u) {
                int p = tree[u].parent;
                int f = tree[p].fail;
                if (tree[f].next[x] != -1) {
                    return tree[f].next[x];
                }
                u = f;
            }
            if (tree[u].next[x] != -1)
                return tree[u].next[x];
            return u;
        };
        while (!q.empty()) {
            int u = q.front();
            q.pop();
            tree[u].fail = compute_failure(u);
            for (int i = 0; i < K; i++) {
                if (tree[u].next[i] != -1) {
                    q.push(tree[u].next[i]);
                }
            }
        }
        for (int i = 0; i < tree.size(); i++) {
            int j = tree[i].fail;
            while (j) {
                if (tree[j].leaf) {
                    break;
                }
                j = tree[j].fail;
            }
            tree[i].exit = j;
        }
    }
    void print() {
        // tree
        cout << "tree:\n";
        for (int i = 0; i < tree.size(); i++) {
            for (int j = 0; j < K; j++) {
                if (tree[i].next[j] != -1) {
                    cout << i << " ";
                    cout << tree[i].next[j] << " ";
                    cout << char(j + 'a') << "\n";

                }
            }
        }
        // failure
        cout << "failure\n";
        for (int i = 0; i < tree.size(); i++) {
            cout << i << " ";
            cout << tree[i].fail << "\n";
        }
        // exit
        cout << "exit:\n";
        for (int i = 0; i < tree.size(); i++) {
            if (tree[i].exit != 0) {
                cout << i << " ";
                cout << tree[i].exit << "\n";
            }
        }
        cout << "\n";
    }

    int fail(int u) { return tree[u].fail; }
    int next(int u, int c) { return tree[u].next[c]; }
    vertex operator[](int u) { return tree[u]; }
};
long long find_strings(const vector<string>& pats, const string& s, vector<bool>& res) {
    // vector<string> words = {
    //     "abba", "baab", "aba", "aabb"
    // };
    // string t = "abbaabaabab";
    // cout << find_strings(words, t) << endl;

    int cnt = 0;
    aho_corasick aho;
    vector<bool> vis(pats.size(), false);

    for (auto& p: pats)
        aho.insert(p, cnt++);
    aho.build_aho();

    int u = 0, ind = -1;
    vector<array<int, 3>> ret; // which string, which index, which node
    for (auto& c : s) {
        ind++;
        int x = c - BASE;
        if (aho.tree[u].next[x] != -1) {
            u = aho.next(u, x);
        }
        else {
            while (u) {
                u = aho.fail(u);
                if (aho.next(u, x) != -1) {
                    break;
                }
            }
            if (aho.next(u, x) != -1) {
                u = aho.next(u, x);
            }
        }

        if (aho[u].leaf) {
            ret.push_back({aho[u].id, ind, u});
            for (auto l : aho[u].ids) {
                vis[l] = true;
            }
        }
        int v = aho[u].exit;
        while (v) {
            ret.push_back({aho[v].id, ind, v});
            for (auto l : aho[v].ids) {
                vis[l] = true;
            }
            v = aho[v].exit;
        }
    }
    res = vis;
    return ret.size();
}

// Clean Code (Yassin Code)
struct aho_corasick {
    int n;
    const int K = 10, BASE = '0';
    vector<int> fail, exit;
    vector<vector<int>> next, out;

    aho_corasick(): n (0) { add_node(); }

    int add_node() {
        fail.emplace_back(0);
        exit.emplace_back(0);
        out.emplace_back();
        next.emplace_back(K, 0);
        return n++;
    }
    int get_idx(char c) { return c - BASE; }
    void insert(const string& s, int idx) {
        int u = 0;
        for (char c: s) {
            if (!next[u][get_idx(c)]) {
                next[u][get_idx(c)] = add_node();
            }
            u = next[u][get_idx(c)];
        }
        out[u].emplace_back(idx);
    }
    void build() {
        queue<int> q;
        q.push(0);

        while (!q.empty()) {
            int u = q.front();
            q.pop();

            for (int i = 0; i < K; i++) {
                int v = next[u][i];
                if (!v) {
                    next[u][i] = next[fail[u]][i];
                }
                else {
                    fail[v] = u ? next[fail[u]][i] : 0;
                    exit[v] = out[fail[v]].empty() ? exit[fail[v]] : fail[v];
                    q.push(v);
                }
            }
        }
    }
    int advance(int u, char c) {
        while (u && !next[u][get_idx(c)])
            u = fail[u];
        u = next[u][get_idx(c)];
        return u;
    }
};

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    return 0;
}
