#include <bits/stdc++.h>
#define int long long
#define all(x) x.begin(), x.end()
#define rall(x) x.rbegin(), x.rend()
#define sz(x) (int)(x).size()

using namespace std;
const double pi = acos(-1);
const int mod = 1e9 + 7;
const int N = 2e5 + 5, M = 1e2 + 2, K = 5e2 + 5;

int dx[] = {0, 0, 1, -1, 1, 1, -1, -1};
int dy[] = {1, -1, 0, 0, 1, -1, 1, -1};

void fast()
{
#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
#endif
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);
}

/*
    Content:
    --- Chapter: Bases ...................... (Done)
    --- Chapter: Polynomials ................ (Done)
    --- Chapter: Sequences .................. (Not Done)
    --- Chapter: Summations ................. (Not Done)
    --- Chapter: Matrices ................... (Not Done)
    --- Chapter: Gaussian Elemenation ....... (Not Done)
    --- Chapter: Polynomial Multiplication .. (Not Done)
*/

/// ========================================================================================================================================================================================================

/// === Chapter: Bases
/*
    --- Base k has k digits: 0 : k-1
    --- Number x in base k = a_n*k^n + a_(n-1)*k^(n-1) + ... + a1*k^1 + a0
    --- There is another representation: 5472 = (((5*10 + 3)*10 + 4)*10 + 7)*10 + 2
    --- x % k = last digits of x in base k
        x / k = removes last digit of x in base k

    --- To convert number in base k to decimal:     --- let: x = (a)(b)(c)(d)...
                                                    --- x = (((a*k + b)*k + c)*k + d)...

    --- To convert number in decimal to base k:     --- let: x = a2*k^2 + a1*k + a0
                                                    --- if we take x%k we get the last digit => a0 = x%k
                                                    --- if we divide by k we remove the last digit => x/k = a2*k + a1
                                                    --- repeate until x != 0
*/
string letters = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
int to_int(char c) { return letters.find(c); }
int from_base_to_decimal(string num, int base)
{
    int res = 0;
    for (int i = 0; i < num.size(); i++)
        res *= 10, res += to_int(num[i]);
    return res;
}
string from_deciaml_to_base(int num, int base)
{
    if (num == 0)
        return "0";
    string res = "";
    while (num)
        res = letters[num % base] + res, num /= base;
    return res;
}

/// ========================================================================================================================================================================================================

/// === Chapter: Polynomials
/*

*/
vector<double> add(vector<double> f1, vector<double> f2)
{
    while (f1.size() < f2.size())
        f1.push_back(0);
    while (f2.size() < f1.size())
        f2.push_back(0);
    vector<double> ans(f1.size(), 0);
    for (int i = 0; i < f1.size(); i++)
        ans[i] = f1[i] + f2[i];
    return ans;
}
vector<double> subtract(vector<double> f1, vector<double> f2)
{
    while (f1.size() < f2.size())
        f1.push_back(0);
    while (f2.size() < f1.size())
        f2.push_back(0);
    vector<double> ans(f1.size(), 0);
    for (int i = 0; i < f1.size(); i++)
        ans[i] = f1[i] - f2[i];
    return ans;
}
vector<double> multiply(vector<double> f1, vector<double> f2)
{
    vector<double> conv(f1.size() + f2.size() - 1, 0);
    for (int i = 0; i < f1.size(); i++)
        for (int j = 0; j < f2.size(); j++)
            conv[i + j] += f1[i] * f2[j];
    return conv;
}
vector<double> division(vector<double> f1, vector<double> f2)
{
    vector<double> deconv(f1.size() - f2.size() + 1, 0);
    int index = f1.size() - f2.size();
    while (f1.size() >= f2.size())
    {
        double subtract = f1.back() / f2.back();
        deconv[index] = subtract;
        index--;
        for (int i = 0; i < f2.size(); i++)
            f1[i + (f1.size() - f2.size())] -= subtract * f2[i];
        while (f1.back() == 0)
            f1.pop_back();
    }
    return deconv;
}
vector<double> reminder(vector<double> f1, vector<double> f2)
{
    while (f1.size() >= f2.size())
    {
        double subtract = f1.back() / f2.back();
        for (int i = 0; i < f2.size(); i++)
            f1[i + (f1.size() - f2.size())] -= subtract * f2[i];
        while (f1.back() == 0)
            f1.pop_back();
    }
    return f1;
}
vector<double> gcd(vector<double> f1, vector<double> f2)
{
    while (f2.size() != 0)
    {
        vector<double> rem = reminder(f1, f2);
        f1 = f2;
        f2 = rem;
    }
    return f1;
}
vector<double> lcm(vector<double> f1, vector<double> f2)
{
    vector<double> lcm = division(multiply(f1, f2), gcd(f1, f2));
    return lcm;
}
vector<double> derivate(vector<double> f)
{
    if (f.size() == 1)
        return vector<double>{0};
    vector<double> ans(f.size() - 1, 0);
    for (int i = 0; i < ans.size(); i++)
        ans[i] = (i + 1) * f[i + 1];
    return ans;
}

/// ========================================================================================================================================================================================================

/// === Chapter: Sequences
/*
    --- Sequence: a list of numbers (in order)
    --- Serie: refer to sequence summation
    --- Sequences can be inreasing, decreasing, increasing-decreasing, positive, negative, or anything
*/

/// ========================================================================================================================================================================================================

int32_t main()
{
    fast();

    return 0;
}