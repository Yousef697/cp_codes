#include <bits/stdc++.h>
#define int long long
#define double long double

using namespace std;

// Sieve
vector<bool> sieve(int n) {
    vector<bool> prime(n + 1, true);
    prime[0] = prime[1] = false;
    for (int i = 4; i <= n; i++) {
        prime[i] = false;
    }
    for (int i = 3; i * i <= n; i += 3) {
        if (prime[i]) {
            for (int j = i * i; j <= n; j += i) {
                prime[j] = false;
            }
        }
    }
    return prime;
}
int count_primes(int n) {
    const int s = 10000;
    vector<int> primes;
    int sq = sqrt(n) + 1;
    vector<bool> isprime(sq + 2, true);
    for (int i = 2; i <= sq; i++) {
        if (isprime[i]) {
            primes.push_back(i);
            for (int j = i * i; j <= sq; j += i) {
                isprime[j] = false;
            }
        }
    }

    int cnt = 0;
    vector<bool> block(s);
    for (int k = 0; k * s <= n; k++) {
        int start = k * s;
        block = vector<bool> (s, true);
        for (auto p : primes) {
            int st = (start + p - 1) / p;
            int j = max(st, p) * p - start;
            for (j; j < s; j += p) {
                block[j] = false;
            }
        }
        if (k == 0) {
            block[0] = block[1] = false;
        }
        for (int i = 0; i < s && start + i <= n; i++) {
            if (block[i]) {
                cnt++;
            }
        }
    }
    return cnt;
}
vector<bool> segmented_sieve(int l, int r) {
    int sq = sqrt(r) + 1;
    vector<int> primes;
    vector<bool> flag(sq + 1, true);
    for (int i = 2; i <= sq; i++) {
        if (flag[i]) {
            primes.push_back(i);
            for (int j = i * i; j <= sq; j += i) {
                flag[j] = false;
            }
        }
    }

    vector<bool> is_prime(r - l + 1, true);
    for (auto i : primes) {
        for (int j = (l + i - 1) / i * i; j <= r; j += i) {
            is_prime[j] = false;
        }
    }
    if (l == 1) {
        is_prime[0] = false;
    }
    return is_prime;
}
vector<int> linear_sieve(int n) {
    vector<int> primes;
    vector<int> spf(n + 1, 0);

    for (int i = 2; i <= n; i++) {
        if (spf[i] == 0) {
            spf[i] = i;
            primes.push_back(i);
        }
        for (auto j : primes) {
            if (j > spf[i] || i * j > n)
                break;
            spf[i * j] = j;
        }
    }
    return spf;
}

// Primarily Tests
int fp(int i, int j, int mod) {
    if (j == 0) {
        return 1;
    }
    int ret = fp(i, j / 2, mod);
    ret = (__int128_t)ret * ret % mod;
    if (j & 1) {
        ret = (__int128_t)ret * i % mod;
    }
    return ret;
}
// Fermat Probabilistic Test
bool prob_fermat(int n, int iter = 10) {
    if (n < 4) {
        return n == 2 || n == 3;
    }
    for (int i = 0; i < iter; i++) {
        int a = 2 + rand() % (n - 3);
        if (fp(a, n - 1, n) != 1) {
            return false;
        }
    }
    return true;
}
// Miller-Rabin Primarily Test
bool check_composite(int n, int a, int d, int s) {
    int x = fp(a, d, n);
    if (x == 1 || x == n - 1) {
        return false;
    }
    for (int i = 1; i < s; i++) {
        x = (__uint128_t)x * x % n;
        if (x == n - 1) {
            return false;
        }
    }
    return true;
}
bool miller_rabin(int n, int iter = 1000) {
    if (n < 4) {
        return n == 2 || n == 3;
    }
 
    int s = 0, d = n - 1;
    while (!(d & 1)) {
        s++;
        d /= 2;
    }
 
    vector<int> primes = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
 
    for (int a : primes) {
        if (n == a) {
            return true;
        }
        if (check_composite(n, a, d, s)) {
            return false;
        }
    }
    return true;
}

// extended gcd & linear diophantine equation
int extended_gcd(int a, int b, int &x, int &y) {
    if (b == 0) {
        x = 1, y = 0;
        return a;
    }
    int x1, y1;
    int g = extended_gcd(b, a % b, x1, y1);
    x = y1;
    y = x1 - y1 * (a / b);
    return g;
}
bool linear_diophantine(int a, int b, int c, int &x, int &y, int &g) {
    g = extended_gcd(abs(a), abs(b), x, y);
    if (c % g)
        return false;

    x *= c / g;
    y *= c / g;
    if (a < 0) x = -x;
    if (b < 0) y = -y;
    return true;
}
void shift_solution(int &x, int &y, int a, int b, int k, int g = 1) {
    x += k * b / g;
    y -= k * a / g;
}
int find_all_solution(int a, int b, int c, int minx, int maxx, int miny, int maxy, int &ans) {
    int x, y, g;
    bool ok = linear_diophantine(a, b, c, x, y, g);
    if (!ok)
        return 0;

    a /= g, b /= g;
    int sign_a = (a > 0 ? 1 : -1);
    int sign_b = (b > 0 ? 1 : -1);

    // find minimum x
    shift_solution(x, y, a, b, (minx - x) / b);
    if (x < minx)
        shift_solution(x, y, a, b, sign_b);
    if (x > maxx)
        return 0;
    int lx1 = x;

    // find maximum x
    shift_solution(x, y, a, b, (maxx - x) / b);
    if (x > maxx)
        shift_solution(x, y, a, b, -sign_b);
    int rx1 = x;

    // find x corresponding to minimum y
    shift_solution(x, y, a, b, -(miny - y) / a);
    if (y < miny)
        shift_solution(x, y, a, b, -sign_a);
    int lx2 = x;

    // find x corresponding to maximum y
    shift_solution(x, y, a, b, -(maxy - y) / a);
    if (y > maxy)
        shift_solution(x, y, a, b, sign_a);
    int rx2 = x;

    if (lx2 > rx2)
        swap(lx2, rx2);

    int lx = max(lx1, lx2);
    int rx = min(rx1, rx2);

    if (lx > rx)
        return 0;
    ans = lx;
    return (rx - lx) / abs(b) + 1;
}

// Chinese Remainder Theorem
int spf[MAXN];
void pre() {
    memset(spf, 0, sizeof(spf));
    for (int i = 1; i < MAXN; i++)
        spf[i] = i;
    for (int i = 2; i < MAXN; i++) {
        if (spf[i] != i)
            continue;

        for (int j = i * i; j < MAXN; j += i)
            if (spf[j] == j)
                spf[j] = i;
    }
}
int inv(int i, int mod) {
    int x, y, g;
    bool ok = linear_diophantine(i, -mod, 1, x, y, g);
    return x;
}
int crt(vector<int> &rems, vector<int> &mods) {
    int M = 1, n = rems.size();
    for (auto m: mods)
        M *= m;

    int ans = 0;
    for (int i = 0; i < n; i++) {
        int N = M / mods[i];
        int N_INV = inv(N, mods[i]);
        ans = (ans + rems[i] * N % M * N_INV + M) % M;
        ans = (ans + M) % M;
    }
    return ans;

    // rems = {4, 3}, mods = {6, 5}, ans = 28 mod 6*5
    // rems = {1, 4, 6}, mods = {3, 5, 7}, ans = 34 mod 3*5*7
}
int general_crt(vector<int> &rems, vector<int> &mods) {
    // run pre function before use this function
    int n = rems.size();
    vector<pair<int, int> > mod_rem;
    for (int i = 0; i < n; i++) {
        int x = mods[i];
        vector<int> p;

        while (x != 1) {
            int y = spf[x], pj = 1;
            while (x % y == 0)
                x /= y, pj *= y;
            p.push_back(pj);
        }
        for (auto pj: p) {
            mod_rem.push_back({pj, rems[i] % pj});
        }
    }

    map<int, pair<int, int> > freq;
    sort(mod_rem.begin(), mod_rem.end());

    for (int i = (int) mod_rem.size() - 1; i >= 0; i--) {
        int x = spf[mod_rem[i].first];
        if (!freq.count(x)) {
            freq[x] = mod_rem[i];
        }
    }

    vector<int> rems2, mods2;
    for (auto &[i, j]: freq) {
        mods2.push_back(j.first);
        rems2.push_back(j.second);
    }

    int ans = crt(rems2, mods2);
    for (int i = 0; i < n; i++) {
        if (ans % mods[i] != rems[i])
            return -1;
    }
    return ans;

    // rems = {3, 5}, mods = {10, 12}, ans = 53 mod lcm(10, 12)
}

// Euler Totient Function
int phis[MAXN];
int phi(int n) {
    int res = n;
    for (int i = 2; i * i <= n; i++) {
        if (n % i)
            continue;

        while (n % i == 0)
            n /= i;
        res -= res / i;
    }
    if (n > 1)
        res -= res / n;
    return res;
}
void all_phi() {
    // n * log(log(n))
    for (int i = 0; i < MAXN; i++)
        phis[i] = i;
    for (int i = 2; i < MAXN; i++) {
        if (phis[i] != i)
            continue;
        for (int j = i; j < MAXN; j += i)
            phis[j] -= phis[j] / i;
    }

    // n * log(n)
    // phis[0] = 0, phis[1] = 1;
    // for (int i = 2; i < MAXN; i++)
    //     phis[i] = i - 1;
    // for (int i = 2; i < MAXN; i++) {
    //     for (int j = 2 * i; j < MAXN; j += i) {
    //         phis[j] -= phis[i];
    //     }
    // }
}

// Mobius Function
int prime[MAXN], mobs[MAXN];
int mobius(int n) {
    int ret = -1;
    for (int i = 2; i * i <= n; i++) {
        int cnt = 0;
        while (n % i == 0)
            n /= i, cnt++;
        if (cnt == 1)
            ret = -ret;
        else if (cnt > 1)
            ret = 0;
    }
    if (n > 1)
        ret = -ret;
    return ret;
}
void all_mobius() {
    for (int i = 0; i < MAXN; i++)
        mobs[i] = -1, prime[i] = 1;
    for (int i = 2; i < MAXN; i++) {
        if (!prime[i])
            continue;
        mobs[i] = 1;
        for (int j = 2 * i; j < MAXN; j += i) {
            prime[j] = 0;
            if (j % (i * i) == 0)
                mobs[j] = 0;
            else
                mobs[j] = -mobs[j];
        }
    }

    // int n = 10, ans = 0, cnt = 0;
    // for (int i = 2; i <= n; i++) {
    //     ans += mobs[i] * (n / i) * (n / i);
    // }
    // for (int i = 1; i <= n; i++) {
    //     for (int j = 1; j <= n; j++) {
    //         if (__gcd(i, j) == 1) {
    //             cnt++;
    //         }
    //     }
    // }
    //
    // cout << ans << "\n";
    // cout << cnt << "\n";
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
    
    int n = 1e6;
    cout << count_primes(n) << "\n";
    auto ret = linear_sieve(20);
    for (int i = 1; i <= 20; i++) {
        cout << i << " " << ret[i] << "\n";
    }
    
    return 0;
}
