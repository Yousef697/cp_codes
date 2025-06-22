// https://vjudge.net/problem/HDU-6222
#include <bits/stdc++.h>

using namespace std;

struct bigint {

    int sign = 1;
    vector<int> digits;

    bigint() {}
    bigint(int v) {
        digits.clear();
        if (v < 0)
            sign = -1, v *= -1;

        while (v)
            digits.push_back(v % 10), v /= 10;
    }
    bigint(string s, int _sign = 1) {
        if (s.size() && s[0] == '-')
            sign = -1, s = s.substr(1);
        sign *= _sign;

        reverse(s.begin(), s.end());
        for (auto i : s)
            digits.push_back(i - '0');
    }
    bigint(const vector<int>& v, int _sign = 1) {
        digits = v;
        sign = _sign;
    }

    int size() const {
        return digits.size();
    }
    int is_zero() const {
        return digits.empty();
    }
    void print() {
        if (size() == 0) {
            cout << 0;
            return;
        }
        if (sign == -1)
            cout << "-";
        for (int i = size() - 1; i >= 0; i--)
            cout << digits[i];
    }

    bigint& operator=(int v) {
        digits.clear();
        while (v)
            digits.push_back(v % 10), v /= 10;
        return *this;
    }
    bigint& operator=(const bigint& other) {
        digits = other.digits;
        return *this;
    }
    bigint& operator=(string s) {
        digits.clear();
        reverse(s.begin(), s.end());
        for (auto i : s)
            digits.push_back(i - '0');
        return *this;
    }

    bool operator<(const bigint& other) const {
        if (sign != other.sign)
            return sign < other.sign;
        if (size() != other.size())
            return (size() < other.size());
        for (int i = size() - 1; i >= 0; i--)
            if (digits[i] != other.digits[i])
                return (digits[i] < other.digits[i]);
        return false;
    }
    bool operator>(const bigint& other) const {
        if (sign != other.sign)
            return sign > other.sign;
        if (size() != other.size())
            return (size() > other.size());
        for (int i = size() - 1; i >= 0; i--)
            if (digits[i] != other.digits[i])
                return (digits[i] > other.digits[i]);
        return false;
    }
    bool operator==(const bigint& other) const {
        if (sign != other.sign)
            return false;
        if (digits.size() != other.digits.size())
            return false;
        for (int i = size() - 1; i >= 0; i--)
            if (digits[i] != other.digits[i])
                return (false);
        return true;
    }

    bigint operator+(const bigint& other) const {
        int n = max(size(), other.size());
        vector<int> res(n);

        if (sign == other.sign) {
            for (int i = 0; i < size(); i++)
                res[i] += digits[i];
            for (int i = 0; i < other.size(); i++)
                res[i] += other.digits[i];
            for (int i = 0; i < n - 1; i++)
                res[i + 1] += res[i] / 10, res[i] %= 10;
            if (res[n - 1] >= 10)
                res.push_back(res[n - 1] / 10), res[n - 1] %= 10;
            while (res.size() && res.back() == 0)
                res.pop_back();
            return bigint(res);
        }
        else if (sign == 1 && other.sign == -1) {
            return *this - other;
        }
        else {
            return other - *this;
        }
    }
    bigint operator+=(const bigint& other) {
        *this = *this + other;
        return *this;
    }

    bigint operator*(const bigint& other) const {
        int n = size();
        int m = other.size();
        vector<int> res(n + m);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < m; j++)
                res[i + j] += digits[i] * other.digits[j];
        for (int i = 0; i < n + m - 1; i++)
            res[i + 1] += res[i] / 10, res[i] %= 10;
        while (res.back() >= 10) {
            int x = res.back();
            res.back() %= 10;
            res.push_back(x / 10);
        }
        while (res.size() && res.back() == 0)
            res.pop_back();
        return bigint(res);
    }
    bigint operator*=(const bigint& other) {
        *this = *this * other;
        return *this;
    }

    bigint operator-(const bigint& other) const {

        int sign = 1;
        int n = max(size(), other.size());
        vector<int> a = digits, b = other.digits;
        while (a.size() < n)
            a.push_back(0);
        while (b.size() < n)
            b.push_back(0);

        if (*this < other)
            sign = -1, swap(a, b);

        for (int i = 0; i < n; i++)
            a[i] -= b[i];

        for (int i = 0; i < n - 1; i++) {
            while (a[i] < 0)
                a[i] += 10, a[i + 1]--;
        }
        while (a.back() < 0)
            a.back() += 10;
        while (a.size() && a.back() == 0)
            a.pop_back();
        return bigint(a, sign);
    }
    bigint operator-=(const bigint& other) {
        *this = *this - other;
        return *this;
    }
};
void print(const bigint& bi) {
    if (bi.is_zero()) {
        cout << 0;
        return;
    }
    if (bi.sign == -1)
        cout << '-';
    for (int i = bi.size() - 1; i >= 0; i--)
        cout << bi.digits[i];
}
void println(const bigint& bi) {
    if (bi.is_zero()) {
        cout << 0;
        return;
    }
    if (bi.sign == -1)
        cout << '-';
    for (int i = bi.size() - 1; i >= 0; i--)
        cout << bi.digits[i];
    cout << "\n";
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    vector<bigint> nums;
    nums.push_back(bigint("4"));
    nums.push_back(bigint("14"));
    nums.push_back(bigint("52"));

    for (int i = 3; i <= 60; i++) {
        bigint res("1");
        res *= nums[i - 1];
        res *= bigint("4");
        res -= nums[i - 2];
        nums.push_back(res);
    }

    int t;
    cin >> t;

    while (t--) {
        string s;
        cin >> s;

        bigint n(s);
        for (int i = 0; i < nums.size(); i++) {
            if (n < nums[i] || nums[i] == n) {
                print(nums[i]);
                break;
            }
        }
    }

    return 0;
}
