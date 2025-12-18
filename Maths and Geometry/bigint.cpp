// https://vjudge.net/problem/HDU-6222
#include <bits/stdc++.h>
#define int long long

using namespace std;

struct bigint {

    int sign = 1;
    vector<int> digits;

    bigint() {
        sign = 1;
        digits.push_back(0);
    }
    bigint(int v) {
        digits.clear();
        if (v < 0)
            sign = -1, v *= -1;
        else
            sign = 1;
        while (v)
            digits.push_back(v % 10), v /= 10;
        normalize(digits);
    }
    bigint(string s, int _sign = 1) {
        if (s.size() && s[0] == '-')
            sign = -1, s = s.substr(1);
        if (s.size() && s[0] == '+')
            sign = +1, s = s.substr(1);
        sign *= _sign;

        reverse(s.begin(), s.end());
        for (auto i : s)
            digits.push_back(i - '0');
        normalize(digits);
    }
    bigint(const vector<int>& v, int _sign = 1) {
        digits = v;
        sign = _sign;
        normalize(digits);
    }

    void normalize(vector<int>& num) const {
        int n = num.size();
        for (int i = 0; i < n - 1; i++) {
            if (num[i] < 0) {
                int b = (abs(num[i]) + 9) / 10;
                num[i] += 10 * b;
                num[i + 1] -= b;
            }
            else {
                num[i + 1] += num[i] / 10;
                num[i] %= 10;
            }
        }
        while (num.size() && num.back() >= 10) {
            int x = num.back() / 10;
            num.back() %= 10;
            num.push_back(x);
        }
        while (num.size() && num.back() == 0)
            num.pop_back();
        if (num.empty())
            num.push_back(0);
    }

    int size() const {
        return digits.size();
    }
    int is_zero() const {
        return size() == 1 && digits[0] == 0;
    }

    bigint& operator=(int v) {
        digits.clear();
        if (v < 0)
            sign = -1, v = -v;
        else
            sign = 1;
        while (v)
            digits.push_back(v % 10), v /= 10;
        return *this;
    }
    bigint& operator=(string s) {
        digits.clear();
        if (s.size() && s[0] == '-')
            sign = -1, s = s.substr(1);
        if (s.size() && s[0] == '+')
            sign = +1, s = s.substr(1);

        reverse(s.begin(), s.end());
        for (auto i : s)
            digits.push_back(i - '0');
        return *this;
    }
    bigint& operator=(const bigint& other) {
        digits = other.digits;
        sign = other.sign;
        return *this;
    }

    bool operator<(const bigint& other) const {
        if (sign != other.sign)
            return sign < other.sign;
        if (size() != other.size())
            return (sign == 1 ? size() < other.size() : size() > other.size());
        for (int i = size() - 1; i >= 0; i--)
            if (digits[i] != other.digits[i])
                return (sign == 1 ? digits[i] < other.digits[i] : digits[i] > other.digits[i]);
        return false;
    }
    bool operator>(const bigint& other) const {
        if (sign != other.sign)
            return sign > other.sign;
        if (size() != other.size())
            return (sign == 1 ? size() > other.size() : size() < other.size());
        for (int i = size() - 1; i >= 0; i--)
            if (digits[i] != other.digits[i])
                return (sign == 1 ? digits[i] > other.digits[i] : digits[i] < other.digits[i]);
        return false;
    }
    bool operator==(const bigint& other) const {
        return !(*this > other || *this < other);
    }
    bool operator<=(const bigint& other) const {
        return *this < other || *this == other;
    }
    bool operator>=(const bigint& other) const {
        return *this > other || *this == other;
    }

    bigint operator-() const {
        bigint res = *this;
        if (!is_zero())
            res.sign *= -1;
        return res;
    }

    bigint operator+(const bigint& other) const {

        int n = max(size(), other.size());
        vector<int> res(n);

        if (sign == other.sign) {
            for (int i = 0; i < size(); i++)
                res[i] += digits[i];
            for (int i = 0; i < other.size(); i++)
                res[i] += other.digits[i];
            normalize(res);
            return bigint(res, sign);
        }
        if (sign == 1 && other.sign == -1) {
            bigint tmp = other;
            tmp.sign *= -1;
            return *this - tmp;
        }
        bigint tmp = *this;
        tmp.sign = -1;
        return other - tmp;
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
        normalize(res);
        return bigint(res, sign * other.sign);
    }
    bigint operator*=(const bigint& other) {
        *this = *this * other;
        return *this;
    }

    bigint operator-(const bigint& other) const {
        int n = max(size(), other.size());
        vector<int> a = digits, b = other.digits;
        while (a.size() < n)
            a.push_back(0);
        while (b.size() < n)
            b.push_back(0);

        if (sign != other.sign) {
            for (int i = 0; i < n; i++)
                a[i] += b[i];
            normalize(a);
            return bigint(a, sign);
        }

        if (*this < other) {
            bigint tmp = other - *this;
            return bigint(tmp.digits, -1);
        }
        for (int i = 0; i < n; i++)
            a[i] -= b[i];
        normalize(a);
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
    cout << endl;
}

int32_t main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    vector<bigint> nums;
    nums.push_back(bigint(4));
    nums.push_back(bigint(14));

    for (int i = 2; i <= 60; i++) {
        int n = nums.size();
        nums.push_back(bigint(4) * nums[n - 1] - nums[n - 2]);
    }

    int t;
    cin >> t;

    while (t--) {
        string s;
        cin >> s;

        bigint n(s);
        for (int i = 0; i < nums.size(); i++) {
            if (n <= nums[i]) {
                println(nums[i]);
                break;
            }
        }
    }

    return 0;
}
