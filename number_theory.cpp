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
    --- Chapter: Modular Arithmetic ............................................... (Done)
    --- Chapter: Primes, Divisors, Factorization .................................. (Done)
    --- Chapter: Factorial ........................................................ (Done)
    --- Chapter: Fibonacci, GCD, LCM, nPr, nCr, Fast Power, Geometric Series Sum .. (Done)
    --- Chapter: Extended GCD ..................................................... (Done)
    --- Chapter: Linear Diophantine ............................................... (Done)
    --- Chapter: Congruence ....................................................... (Done)
    --- Chapter: Euler Totient Function ........................................... (Done)
    --- Chapter: Mobius Function .................................................. (Done)
    --- Chapter: Modular Inverse .................................................. (Done)
    --- Chapter: Modular Arithmetic Applications .................................. (Done)
    --- Chpater: Chinese Remainder Theorem ........................................ (Done)
    --- Chpater: Permutations ..................................................... (Done)
*/

/// ========================================================================================================================================================================================================

/// === Chpater: Modular Arithmetic
/*
    --- a modulo n: finds the reminder of a / n [~ a % n]
        a can be written as: x*n + r, 0 <= |r| < n => r is a % n
        r = a % n = a - floor(a / n) * n
        0 < |r| < n ( |r| has n values: [0, n-1] )

    --- -ve is not the same as +ve
        -10 % 7 is not equal to 10 % 7
        r = a % n = (a + q*n) % n
        so -10 % 7 = (-10+7) % 7 = -3 % 7 = -3
        taking the mod and adding 1 cycle and taking the mod always gives a positive answer
        (a%n + n) % n is always positive

    --- mod and division are time expensive operations, so we can write an iterative code
        while (a >= n) a -= n;
        while (a < 0) a += n;
        (we add or remove cycles to get modulo)

    --- a % n = 0  =>  a is divisible by n
        a % n = b % n  =>  (a - b) % n = 0
        largest n such that a % n = b % n  =>  n = b - a

    --- Some Facts: --- (a % n) % n = a % n
                    --- (n ^ x) % n = 0 [x > 0, x is integer]
                    --- -a % n != a % n
                    --- (a + b) % n = (a%n + b%n) % n
                    --- x % (a + b) != x%a + x%b
                    --- x%10: last digit, floor(x/10): removes last digit
                    --- (a ^ b) % n = ( (a % n) ^ b ) % n
                    --- (1 / a) % n  =>  Mod multiplicative inverse
                    --- ( (a*b)%n * (1/a)%n ) % n = b % n
                    --- a % (2 ^ n) = a & (2^n-1)
*/
int division_mod(int a, int b) { return a - a / b * b; }
int positive_mod(int a, int b) { return (a % b + b) % b; }
int iterative_mod(int a, int b)
{
    while (a >= b)
        a -= b;
    while (a < 0)
        a += b;
    return a;
}
int add_mod(int a, int b, int mod) { return (a % mod + b % mod) % mod; }
int last_digit(int a) { return a % 10; }
int power_mod(int a, int b, int mod) { return (int)pow(a % mod, b) % mod; }
int power_2_mod(int a, int mod) { return (a & (mod - 1)); }

/// ========================================================================================================================================================================================================

/// === Chpater: Primes, Divisors, Factorization
/*
    --- Prime Checking: --- Main Logic: Check all numbers from 2:n-1
                        --- Optimization 1: check if the number is 2 or its multiples
                                            then check odd numbers
                        --- Optimization 2: check numbers from 2:sqrt(n)
                                            becuase if n = a * b the a <= sqrt(n), b >= sqrt(n)
                        --- Optimization 3: calculate sqrt(n) before enterig the loop

    --- Prime Seive:    --- for each number i from 2 to n: if i was not marked, marke all multiples of i

    --- Factorization:  Getting all divisors of a number

    --- Generating Divisors:    --- loop on numbers i from 1 to sqrt(n)
                                --- if n mod i = 0, add i and n / i to the prime factors
                                --- after the loop check if n has an integer square root or not

    --- Prime Factorization:    --- loop on numbers i from 1 to sqrt(n)
                                --- while n mod i = 0, add i to the divisors and then divide n by i
                                --- after the loop check if n is not 1, add n to the prime factors

    --- Number of Divisors: --- assume n = p1^a1 * p2^a1 * ... * pk^ak
                            --- number of divisors is (a1 + 1) * (a2 + 1) * ... * (ak + 1)

    --- Number of divisors of n^z:  --- assume n = p1^a1 * p2^a1 * ... * pk^ak
                                    --- n^z = p1^(a1*z) * p2^(a1*z) * ... * pk^(ak*z)
                                    --- number of divisors is (a1*z + 1) * (a2*z + 1) * ... * (ak*z + 1)

    --- Number of diviors using Seive:  --- for number i, incearse number of divisors of all multiples if i by 1
*/
vector<int> primes(N, 1);
bool is_prime(int n)
{
    if (n <= 1)
        return false;
    if (n == 2 || n == 3)
        return true;
    for (int i = 2; i * i <= n; i++)
    {
        if (n % i == 0)
            return false;
    }
    return true;
}
void all_prime()
{
    // memset(all(primes), 1, sizeof primes);
    for (int i = 0; i < N; i++)
        primes[i] = 1;

    primes[0] = 0, primes[1] = 0;
    for (int i = 2; i * i < N; i++)
    {
        if (!primes[i])
            continue;

        for (int j = i * i; j < N; j += i)
            primes[j] = 0;
    }
}
vector<int> divisors(int n)
{
    vector<int> ans;
    int i;
    for (i = 1; i * i < n; i++)
        if (n % i == 0)
            ans.push_back(i), ans.push_back(n / i);

    if (i * i == n)
        ans.push_back(i);
    return ans;
}
vector<int> prime_factorization(int n)
{
    vector<int> ans;
    for (int i = 2; i * i <= n; i++)
        while (n % i == 0)
            ans.push_back(i), n /= i;

    if (n != 1)
        ans.push_back(n);

    return ans;
}
vector<pair<int, int>> prime_factorization_power_form(int n)
{
    vector<pair<int, int>> ans;
    for (int i = 2; i * i <= n; i++)
    {
        int cnt = 0;
        while (n % i == 0)
            cnt++, n /= i;

        if (cnt != 0)
            ans.push_back({i, cnt});
    }
    if (n != 1)
        ans.push_back({n, 1});

    return ans;
}

/// ========================================================================================================================================================================================================

/// === Chapter: Factorial
/*
    --- Factorial of x: result of multiplication of all numbers between 1 and x

    --- Factorial x is x!

    --- To get power of x inside n!: let m = n!, while m%x=0: (m /= x, number_of_powers++)

    --- Some facts: --- 0! = 1! = 1
                    --- n! % x = 0 for 1 <= x <= n
                    --- (p-1)! % p = p-1 iff p is prime (Wilson's Theorem)
                    --- Number of digits of n! = 1 + floor( sum of log_10(i) ) for 1 <= i <= n;
                    --- Number of trailing zeros: number of power of 5 inside n!
                    --- Right most non-zero digit is:   --- (n! / 10^number of trailing zeros)%10
                                                        --- or delete 10's from the prime factorization of n!

    --- Given m, find smallest n such that n! has m trailing zeros
*/
vector<int> fact(21);
int factorial(int x)
{
    if (x == 0 || x == 1)
        return 1;
    int ans = 1;
    for (int i = 1; i <= x; i++)
        ans *= i;
    return ans;
}
int recursive_factorial(int x)
{
    if (x == 0 || x == 1)
        return 1;
    return x * recursive_factorial(x - 1);
}
int number_of_digits_of_factorial_x(int x)
{
    long double ans = 0;

    for (int i = 1; i <= x; i++)
        ans += log10(i);
    return (int)ans + 1;
}
int power_of_x_in_factorial_n(int x, int n)
{
    int a = x, ans = 0;

    while (a <= n)
        ans += n / a, a *= x;
    return ans;
}
int number_of_trailing_zeros_in_factorial_n(int n)
{
    return power_of_x_in_factorial_n(5, n);
}
int right_most_non_zero_digit_in_factorial_n(int n)
{

    int ans = 1;

    int x = power_of_x_in_factorial_n(5, n);
    int y = power_of_x_in_factorial_n(2, n) - x;

    while (y--)
        ans *= 2, ans %= 10;

    for (int i = 3; i <= n; i++)
    {
        if (i == 5)
            continue;
        if (!primes[i])
            continue;

        int a = power_of_x_in_factorial_n(i, n);

        while (a--)
            ans *= i, ans %= 10;
    }
    return ans;
}
void all_fact()
{
    fact[0] = fact[1] = 1;
    for (int i = 2; i <= 10; i++)
        fact[i] = i * fact[i - 1];
}
/// ========================================================================================================================================================================================================

/// === Chapter: Fibonacci, GCD, LCM, nPr, nCr, Fast Power, Geometric Series Sum, Stirling Numbers
/*
    --- Fibonacci:  --- Each number is the sum of the two previous numbers
                    --- The base cases: fib(0) = 0, fib(1) = 1
                    --- fib(n) = fib(n-1) + fib(n-2)

    --- GCD:    --- Greatest common divisor that divides two numbers a, b
                --- let n = gcd(a, b):  --- a%n = b%n = 0
                                        --- (a+b)%n = (a-b)%n = 0
                                        --- (a%n + b%n)%n = (0 + 0)%n = 0
                                        --- gcd(a, b) = gcd(a-b, b) = gcd(a-2b, b) = ...
                                        --- gcd(a, b) = gcd(a%b, b) = gcd(b, a%b);
                --- If you have the prime factorization of a, b, the gcd is
                    the minimum power of p for all common primes inside a, b

    --- LCM:    --- Least common multiple that a, b divide it
                --- If you have the prime factorization of a, b, the lcm is
                    the maximum power of p for all primes inside a, b
                --- lcm(a, b) = a / gcd(a, b) * b

    --- Permutation:    --- Number of arrangements r objects out of n objects (we care about order)
                        --- It depends on the rule of product:
                            First you have n choices then n-1 the n-2 then ... then n-r+1
                        --- P(n, r) = n! / (n-r)!

    --- Combinations:   --- Number of arrangements r objects out of n objects (we don't care about order)
                        --- C(n, r) = n! / ( (n-r)! * r! ) = P(n, r) / r!

    --- Fast Power: --- If you have 5^6, you can calculate 5^3 and square it
                    --- If you have 5^9, you can calculate 5^4 and square it and multiply by 5

    --- Geometric Series Sum:   --- If you have a^1 + a^2 + a^3 + a^4 + a^5 + a^6
                                --- You can simplify: a^1 + a^2 + a^3 + a^3(a^1 + a^2 + a^3)
                                --- Again: (a^1 + a^2 + a^3) * (1 + a^3)
                                --- Again: (a^1 + a^2 + a^3) * (1 + a^3 + a^2 + a^1 - (a^2 + a^1))
                                --- Generalization: if k even: ans(k) = ans(k/2) * (1 + ans(k/2) - ans(k/2 - 1))
                                                    if k odd : ans(k) = a * (1 + ans(k-1))
                                                    if k = 0 : return 0
*/
int fibonacci(int n)
{
    if (n <= 1)
        return n;
    int a = 0, b = 1, c;
    for (int i = 2; i <= n; i++)
        c = a + b, a = b, b = c;
    return b;
}
int gcd(int a, int b)
{
    int temp;
    while (b)
        temp = a, a = b, b = temp % a;
    return a;
}
int recursive_gcd(int a, int b)
{
    if (b == 0)
        return a;
    return gcd(b, a % b);
}
int lcm(int a, int b)
{
    return a / gcd(a, b) * b;
}
int permutaions(int n, int k)
{
    return factorial(n) / factorial(n - k);
}
int combinations(int n, int k)
{
    return factorial(n) / (factorial(n - k) * factorial(k));
}
int fast_power(int n, int p)
{
    if (p == 0)
        return 1;

    int ans = fast_power(n, p / 2);
    ans *= ans;

    if (p % 2 == 0)
        return ans;
    else
        return ans * n;
}
int fast_power_mod(int n, int p, int mod)
{
    if (p == 0)
        return 1;

    int ans = fast_power_mod(n, p / 2, mod);
    ans = (ans % mod) * (ans % mod);

    if (p % 2 == 0)
        return ans % mod;
    else
        return (ans % mod * n % mod) % mod;
}
int geometric_series_sum(int a, int n)
{
    if (n == 0)
        return 0;

    if (n % 2 == 1)
        return a * (1 + geometric_series_sum(a, n - 1));

    int x = geometric_series_sum(a, n / 2);
    int y = geometric_series_sum(a, n / 2 - 1);

    return x * (1 + x - y);
}

/// ========================================================================================================================================================================================================

/// === Chapter: Extended GCD
/*
    --- Extended GCD solve equations like: x*a + y*b = gcd(a, b) for x, y

    --- The code is like GCD code except for:
        --- if b = 0, we can solve the eqaution by putting x=1, y=0 so: 1*a + 0*0 = a
        --- so we now solved the equation for gcd(b, a%b) now we want to solve it for gcd(a, b)
        --- remember:   gcd(a, b) = x*a + y*b
                        gcd(b, a%b) = x1*b + y1*(a%b) = gcd(a, b) [a%b = a - floor(a/b) * b]
                        gcd(a, b) = x1*b + y1*(a - floor(a/b) * b)
                        gcd(a, b) = y1*a + (x1 - floor(a/b) * y1)*b
        --- Then:   The next x = the previous y
                    The next y = the previous x - floor(a/b) * the previous y
        --- We can generate other solutions:
                    (xo, yo) = (x + k*b/gcd(a, b), y - k*a/gcd(a, b)) for all values of k
*/
vector<int> extended_gcd(int a, int b)
{
    if (b == 0)
        return vector<int>{a, 1, 0};

    vector<int> v = extended_gcd(b, a % b);

    int x = v[2], y = v[1] - v[2] * a / b;

    return vector<int>{v[0], x, y};
    // Other solutions: (xo, yo) = (x + k*b/gcd(a, b), y - k*a/gcd(a, b)) for all values of k
}

/// ========================================================================================================================================================================================================

/// === Chapter: Linear Diophantine
/*
    --- Diophantine equations:  --- WHO CARES????
                                --- Just jokking, use extended GCD to solve: a*x + b*y = c
                                --- If c = gcd(a, b), the answer from the extended GCD is true
                                --- else solve for c = gcd(a, b) then multiply by c / gcd(a, b)
                                --- The algorithm work iff c % gcd(a, b) = 0
*/
pair<int, int> linear_diophantine(int a, int b, int c)
{
    vector<int> ans = extended_gcd(a, b);

    if (c % ans[0] != 0)
        return {-1e9, -1e9};

    return {ans[1] * c / ans[0], ans[2] * c / ans[0]};
}

/// ========================================================================================================================================================================================================

/// === Chapter: Congruence
/*
    --- a and b are called congruent modulo n (written a = b (mod n))
        --- It means that: a%n = b%n = x;
        --- (a-b)%n = 0, a-b = q*n

    --- if (a*x = a*y (mod n)) and gcd(a, n) = d, we can rearrange the equation to (x = y (mod n/d))

    --- Some Facts: --- a = b (mod n), then a^m = b^m (mod n), m >= 1
                    --- (x + y)^p = x^p + y^p (mod p), p is prime
                    --- a = b (mod n) and b = c (mod n) then a = c (mod n)
                        a = b (mod n) and c = d (mod n) then a+c = b+d (mod n)
                        a = b (mod n) and c = d (mod n) then a*c = b*d (mod n)
                        a = b (mod n) then a+c = b+c (mod n)
                        a = b (mod n) then a*c = b*c (mod n)
                    --- a*x = b (mod n) then x = b*a^-1

    --- Large Powers:   --- assume that we want to calculte: a^m mod n
                        --- we want to reduce the large power
                        --- we want to find such k that a^k mod n = -1, 0, 1 (small values)
                        --- rewrite the equation: a^(q*k + r) mod n
                        --- ((a^k)^q * a^r) mod n

    --- Linear Modular Equation:    --- solve: a*x = b (mod m)
                                    --- rewrite the equation: a*x - b = y*m
                                    --- ax + y*m = b (mod m)
                                    --- Note that: all possible solution ate between 0, m-1 (due to properties of modular arithmetic)
*/
vector<int> linear_modular_equation(int a, int b, int mod)
{

    vector<int> v = extended_gcd(a, mod), sols;
    int g = v[0], x = v[1];

    if (b % g != 0)
        return sols;

    x = ((x * b / g) % mod + mod) % mod; // In case of gcd(a, mod) != 1

    for (int i = 0; i < g; i++) // Bezout Identity
        sols.emplace_back((x + i * mod / g) % mod);
    sort(all(sols));
    return sols;
}

/// ========================================================================================================================================================================================================

/// === Chapter: Euler Totient Function
/*
    ---Euler Totient Function : --- Phi Function : counts the number of coprimes to n
                                --- Comprimes : gcd(a, b) = 1
                                --- phi(a * b * c) = phi(a) * phi(b) * phi(c) iff a, b, c, are coprimes
                                --- phi(p^k) = p^k - p^(k-1) = p^(k-1) * (p-1)
                                --- phi(1) = phi(2) = 1
                                --- phi(n) is even for n > 2
                                --- phi(n^k) = n^(k-1) * phi(n)
                                --- sqrt(n) <= phi(n) <= n - sqrt(n) except for 2, 6
                                --- n = sum of phi(i) for all i divides n
                                --- phi(n) * d(n) = sum of gcd(k-1, n) for all k such that gcd(k, n) = 1
                                    and d(n) is number of divisors of n
                                --- We can brute force to calculate phi(n) or we can do better
                                --- We can represent n as: n = p1^a1 * p2^a2 * ... * pk^ak
                                --- p1, p2, ..., pk are coprimes
                                --- so: phi(n) = phi(p1^a1) * phi(p2^a2) * ... * phi(pk^ak)
                                --- phi(n!) = product of (if i prime ? i-1 else i) for all 2 <= i <= n
                                --- Totient Range Query:    --- loop on all prime numbers i from 2 to n
                                                            --- make totient[i] = n-1
                                                            --- loop on all multiples of i in the range
                                                            --- multiply totient[j*i] by the answer
*/
vector<int> phis(N);
int phi_of_prime_power_k(int p, int k)
{
    int ans = fast_power(p, k);
    return ans - ans / p;
}
int phi(int n)
{
    if (n == 1 || n == 2)
        return 1;

    vector<pair<int, int>> v = prime_factorization_power_form(n);
    int ans = 1;
    for (auto &[i, j] : v)
        ans *= phi_of_prime_power_k(i, j);
    return ans;
}
int phi_of_n_power_k(int n, int k)
{
    return fast_power(n, k - 1) * phi(n);
}
int phi_of_factorial_n(int n)
{
    int ans = 1;
    for (int i = 2; i <= n; i++)
    {
        ans *= (primes[i] ? i - 1 : i);
    }
    return ans;
}
void all_phi()
{
    // memset(all(phis), 1);
    for (int i = 0; i < N; i++)
        phis[i] = 1;

    for (int i = 2; i < N; i++)
    {
        if (!primes[i])
            continue;

        phis[i] = i - 1;

        for (int j = 2 * i; j < N; j += i)
        {
            int cnt = 0, x = j;
            while (x % i == 0)
                cnt++, x /= i;
            phis[j] *= phi_of_prime_power_k(i, cnt);
        }
    }
}

/// ========================================================================================================================================================================================================

/// === Chapter: Mobius Function
/*
    --- Square-free Integer:    --- A number which there is no repeated prime factor in it (all powers <= 1)
    --- Mobius Function:    --- mobius(n) [μ(n)] =  [
                                                        0: n is not square free
                                                        +1: n is one or n is square free and has even number of prime factors
                                                        -1: n is square free and has odd number of prime factors
                                                    ]
                            --- Mainly, it is used in inclusion-exclusion problems, like this
                            --- It can be used to calculate the index of a number in square-free list
                                let x be the value, let y be the index and it is initially is x
                                for all numbers i from 2 to sqrt(x)
                                y := y + mobius[i] * ( x / (i*i) )
                                [
                                    --- if i = 2 it will remove 4, 8, 12, 16, ..
                                    --- if i = 3 it will remove 9, 18, 27, 81, ..
                                    --- if i = 4 it will not remove anything because 2 has done it before
                                    --- if i = 5 it will remove 25, 50, 75, 100, ..
                                    --- if i = 6 it will remove 36, 72, 144, 288, ...
                                            wait, 36, 72, 144, 288, ... are removed twice by 2 and 3
                                            so they must be add because they are deleted twice
                                ]
                            --- Count number of triple such that: a, b, c <= n and gcd(a, b, c) = 1
                                --- Lets do a reverse thinking: the answer is n^3 - {number of triple such that: a, b, c <= n and gcd(a, b, c) != 1}
                                --- let the answer be n^3
                                --- for all numbers i from 2 to n
                                --- answer := answer + mobius[i] * (n/i)^3 [Why?]
                                [
                                    --- if i = 2 there are n/2 multiple of 2 in the range so for all of these number
                                        the gcd is not 1, we will remove (n/2)^3 from the answer
                                    --- if i = 2 there are n/3 multiple of 3 in the range so for all of these number
                                        the gcd is not 1, we will remove (n/3)^3 from the answer
                                    --- if i = 4 the answer will not change because 2 has done it before
                                    --- if i = 5 there are n/5 multiple of 5 in the range so for all of these number
                                        the gcd is not 1, we will remove (n/5)^3 from the answer
                                    --- if i = 6 there are n/6 multiple of 6 in the range so for all of these number
                                        the gcd is not 1, we will remove (n/6)^3 from the answer
                                            wait, the answer of 6 has been removed twice by 2 and 3
                                            so the answer of 6 must be added to the answer
                                ]
                            --- Mobius Range Query: --- initialize mobius[i] = 1
                                                    --- loop on all prime numbers i from 2 to n
                                                    --- mobius[i] = 1
                                                    --- loop on all multiples of i in the range
                                                    --- mobius[j] = (j % (i*i) == 0 ? 0 : -mobius[j])
*/
vector<int> mobiuses(N);
pair<bool, int> square_free(int n)
{
    vector<pair<int, int>> v = prime_factorization_power_form(n);

    bool ok = 1;
    for (auto &[i, j] : v)
        if (j > 1)
            ok = 0;

    return {ok, v.size()};
}
int mobius(int n)
{
    pair<bool, int> p = square_free(n);
    if (!p.first)
        return 0;
    return (p.second % 2 == 0 ? 1 : -1);
}
void all_mobius()
{
    // memset(all(mobiuses), 1);
    for (int i = 0; i < N; i++)
        mobiuses[i] = 1;

    for (int i = 2; i < N; i++)
    {
        if (!primes[i])
            continue;

        mobiuses[i] = 1;
        for (int j = 2 * i; j < N; j += i)
            mobiuses[j] = (j % (i * i) == 0 ? 0 : -mobiuses[j]);
    }
}

/// ========================================================================================================================================================================================================

/// === Chapter: Modular Inverse
/*
    --- Modular Multiplicative Inverse?! Modular Multiplicative Inverse?? Modular Multiplicative Inverse!!
    --- Modular Multiplicative Inverse: a * (what) = 1 (mod m)
    --- We can calculate (what) iff gcd(a, m) = 1
    --- To find Modular Multiplicative Inverse of a considering m, you want to find a^-1 mod m
    --- Follow my leads:    --- a*x = 1 (mod m)
                            --- a*x - 1 = y*m
                            --- a*x + m*(-y) = 1 (Wait, Extended GCD? YES!!)
                            --- The x is the answer (a^-1);
    --- Follow my leads again:  --- Euler has a theorem: if gcd(a, m) = 1, a^phi(m) = 1 (mod m)
                                --- a^(phi(m) - 1) = 1/a = a^-1 (mod m)
                                --- if m is prime, a^-1 = a^(m-2) (mod m)
    --- Mod Inverse for range:  --- Given p, calculate mod inverse of all numbers 1 from 1 to p-1
                                --- Follow my leads again:  --- p%i = p - (p/i)*i
                                                                since i < p, (p%i)%p = p%i
                                                            --- (p%i)%p = p%p - [(p/i) * i]%p
                                                            --- p%i = - (p/i) * i (mod p)
                                                                divide by (i * p%i)
                                                            --- 1/i = - (p/i) * 1/(p%i) (mod p)
                                                            --- inv[i] = - (p/i) * inv[p % i] (mod p)
                                                            --- inv[i] = p - (p/i) * inv[p % i] (mod p)
                                                                [initailize all values of inv by 1]
    --- Euler Theorem and LARGE Powers: --- remember:   a^phi(m) = 1 (mod m), a^(phi(m)-1) = 1/a (mod m) [gcd(a, m) = 1]
                                                        a^(p-1) = 1 (mod p), a^(p-2) = 1/a (mod p) [p is prime]
                                        --- a^b = a^(n*phi(m) + k)  = a^(n*phi(m)) * a^k (mod m)
                                                                    = [a^phi(m)]^n * a^k (mod m)
                                                                    = [1]^n * a^k (mod m)
                                                                    = a^k (mod m)
                                                                    = a^(b % phi(m)) (mod m)
                                        --- Compute (1/a^b) mod p   = (1/a % p)^b mod p
                                                                    = (a^(p-2) % p)^b mod p
                                                                        [Use the above relation]
                                                                    = a^[(p-2)%(p-1) * b%(p-1)] mod p
                                                                    = a^[-1 * b%(p-1)] mod p
                                                                    = a^[p-1 - * b%(p-1)] mod p
*/
vector<int> invs(N);
int mod_inv(int a, int n)
{
    vector<int> v = extended_gcd(a, n);

    return v[1];
}
int mod_inv_euler(int a, int n)
{
    return fast_power_mod(a, phi(n) - 1, n);
}
int mod_inv_prime(int a, int n)
{
    return fast_power_mod(a, n - 2, n);
}
void all_mod_inv(int p)
{
    for (int i = 0; i < N; i++)
        invs[i] = 1;
    for (int i = 2; i < N; i++)
        invs[i] = p - p / i * invs[p % i], invs[i] = (invs[i] % p + p) % p;
}
int euler_theorem_for_large_powers(int a, int b, int mod)
{
    return fast_power_mod(a, b % phi(mod), mod);
}

/// ========================================================================================================================================================================================================

/// === Chapter: Modular Arithmetic Applications
/*
    --- Factoial mod p after removing all p in the factorial:
                    --- Assume we want to calculate 33! % 5 after removing all p in the factorial
                    --- 1 * 2 * 3 * 4 * 5 * 6 * 7 * 8 * 9 * 10 * 11 * ...
                        1 * 2 * 3 * 4 * 1 * 6 * 7 * 8 * 9 *  2 * 11 * ...
                        So the problem will be
                    ---  1  2  3  4  1 ( 5)
                            6  7  8  9  2 (10)
                        11 12 13 14  3 (15)
                        16 17 18 19  4 (20)
                        21 22 23 24  5 (25)
                        26 27 28 29  6 (30)
                        31 32 33
                    --- Then we take mode 5 for all numbers
                            1  2  3  4  1 ( 5)
                            1  2  3  4  2 (10)
                            1  2  3  4  3 (15)
                            1  2  3  4  4 (20)
                            1  2  3  4  5 (25)
                            1  2  3  4  6 (30)
                            1  2  3
                    --- We now har this: (4! % 5)^6 * (3! % 5) * (a smaller problem with 6 instead of 33 F(6))
                    --- Recall: (p-1)! % p = p-1 = -1 (Wilson's Theorem)
                    --- Then: (-1)^6 * (3! % 5) * (F(6) % 5)

    --- Combinations mod p:     --- n choose k = n! / (k! * (n-k)!)
                                --- let: a = power of p in n!, b = power of p in k!, c = power of p in (n-k)!
                                --- if a-b-c != 0 then the combinations is divisible by p
                                --- remove all powers of p in n!, k!, and (n-k)!
                                --- then for each factorial calculate mod p using the above algorithm after removing all the powers

    --- Locus Theorem:  --- to calculate n choose k mod p
                        --- first write n and k in base p
                        --- let: n = n_k*p^k + n_(k-1)*p^(k-1) + ... + n1*p^1 + n0
                                 m = m_k*p^k + m_(k-1)*p^(k-1) + ... + m1*p^1 + m0
                        --- (n choose k) mod p = (n_k choose m_k) * (n_(k-1) choose m_(k-1)) * ... * (n1 choose m1) * (n0 choose m0)
                        --- There is a theorem called Generalized Locus Theorem that handles mod p^k
                        --- If the mod is not prime, we can write the mod as: p1^k1 * p2^k2 * ...
                        --- calculate the combination using the generalized locus theorem for p_i^k_i
                        --- Using the Chinese Remainder Theorem, We can solve the main question

    --- Catalan Numbers:    --- Cn = (2n choose n) / (n+1) = (2n)! / (n! * (n+1)!) = (2n choose n) - (2n choose (n+1))
                            --- It is all about combinations and factorials, we can solve it using the above algorithms
*/
int factorial_mod_p(int n, int p)
{
    int res = 1;
    while (n)
    {
        for (int i = 1; i <= n % p; i++)
            res = (res * i) % p;
        n /= p;
        if (n % 2 == 1)
            res = p - res;
    }
    return res;
}
int combinations_mod_p(int n, int k, int p)
{
    int a = power_of_x_in_factorial_n(p, n);
    int b = power_of_x_in_factorial_n(p, k);
    int c = power_of_x_in_factorial_n(p, n - k);

    if (a - b - c)
        return 0;

    int up = factorial_mod_p(n, p);
    int down = (factorial_mod_p(k, p) * factorial_mod_p(n - k, p)) % p;
    int down_inv = mod_inv_prime(down, p);

    return (up * down_inv) % p;
}
int combinations_mod_p_locus(int n, int k, int p)
{
    auto to_base_p = [&](int x, int p)
    {
        vector<int> num;
        while (x)
            num.emplace_back(x % p), x /= p;
        return num;
    };

    vector<int> np = to_base_p(n, p), kp = to_base_p(k, p);
    while (kp.size() < np.size())
        kp.emplace_back(0);
    int res = 1;
    for (int i = 0; i < np.size(); i++)
    {
        for (int j = 0; j < k; j++)
            (res *= n - j) %= p;

        int a = factorial_mod_p(k, p);
        a = mod_inv_prime(a, p);

        res = (res % p * a % p) % p;
    }
    return res;
}
int catalan(int n, int p)
{
    int a = combinations_mod_p(2 * n, n, p);
    int b = combinations_mod_p(2 * n, n + 1, p);
    return (a - n + p) % p;
}

/// ========================================================================================================================================================================================================

/// === Chapter: Chinese Remainder Theorem
/*
    #11: Chinese Remainder Theorem
    --- This theorem solves system of equations like this:  x = a1 (mod n1)
                                                            x = a2 (mod n2)
                                                            ...
                                                            x = ak (mod nk)
    --- Brute Force? Bad. Lets get into the theorem
    --- To solve the system, there are some conditions must holds:  --- n1, n1, ..., nk must be pairwise coprime
                                                                    --- ai = aj (mod gcd(ni, nj)) for all i, j
                                                                    --- x are then congreunt to lcm(n1, n2, ..., nk)
    --- if a, b, c are pairwise coprime, then gcd(a, b) = gcd(b, c) = gcd(a, c) = 1, and lcm(a, b, c) = a * b * c
    --- To covert a%n to b%n, just make it as (a/a * b)%n = (a * a^-1 * b)%n
    --- Now the system is: a[i] and mods[i] and b[i] and c[i] and answer=0
    --- Follow my leads:    --- Step 1: for the ith equation, calculate the product of all mods except mods[i]
                                and then assign it to b[i] then take mod mods[i] and assign it to c[i]
                            --- now we want to convert b[i] to c[i] (mod mods[i])
                            --- using the relationship above: (b[i] * b[i]^-1 * c[i]) (mod mods[i]) and add it to the answer
                            --- Boom, we solved it, x is the answer [we can take mod lcm(mods) to get the smallest answer]

    --- Now we want to solve if the mods are not pairwise coprimes
    --- if we have two equations, can we solve it? lets try
    --- Follow my leads:    --- The equations are: x = a1 mod n1, x = a2 mod n2
                            --- x = b1*n1 + a1 = b2*n2 + a2
                            --- n1*b1 + n2*(-b2)x = a2-a1 (wiat, Diophantine? YES!!);
                            --- new rem: rem + mod*x, new mod: lcm(mod, mods[i]), rem %= mod

    --- There is another way to solve if mods are not pairwise coprime
    --- if n1, n2 are not coprime? extend the 2 equetions to 4 equations
        [let d = gcd(n1, n2)]
        x = (a1 % d) mod d
        x = (a1 % (n1/d)) mod n1/d
        x = (a2 % d) mod d
        x = (a2 % (n2/d)) mod n2/d

    --- if we want to compute some value x % c where c is not prime?
    --- Write c as p1^a1 * p2^a2 * ... * pk^ak = b1 * b2 * ... * bk
    --- b1, b2, ..., bk are pairwise coprimes
    --- Now compute x mod each value of b
    --- We want x % c? Chinese Remainder Theorem can solve this system

    --- If we want to calc f(), but there is an intermediate overflow
    --- let: M = 257 * 263 * 269 * 271 = p1 * p2 * p3 * p4
    --- Calculate f() % pi and then using crt we can get the value of f()

*/
int crt_coprime_mods(vector<int> &rems, vector<int> &mods)
{
    int product = 1, ans = 0;
    for (auto &i : mods)
        product *= i;

    for (int i = 0; i < rems.size(); i++)
    {
        int sub = product / mods[i];
        ans += sub * mod_inv_euler(sub, mods[i]) * rems[i];
    }
    return ans % product;
}
int general_crt(vector<int> &rems, vector<int> &mods)
{
    int rem = rems[0], mod = mods[0];

    for (int i = 1; i < rems.size(); i++)
    {
        int a = mod, b = -mods[i], c = rems[i] - rem;
        auto [x, y] = linear_diophantine(a, b, c);

        if (x == -1e9 && y == -1e9)
            return -1;

        rem += mod * x;
        mod = lcm(mod, mods[i]);
        rem = (rem % mod + mod) % mod;
    }
    return rem;
}
int composite_mod(int a, int mod) {}

/// ========================================================================================================================================================================================================

/// === Chapter: Permutations
/*
    index 00: 0 1 2 3
    index 01: 0 1 3 2
    index 02: 0 2 1 3
    index 03: 0 2 3 1
    index 04: 0 3 1 2
    index 05: 0 3 2 1

    index 06: 1 0 2 3
    index 07: 1 0 3 2
    index 08: 1 2 0 3
    index 09: 1 2 3 0
    index 10: 1 3 0 2
    index 11: 1 3 2 0

    index 12: 2 0 1 3
    index 13: 2 0 3 1
    index 14: 2 1 0 3
    index 15: 2 1 3 0
    index 16: 2 3 0 1
    index 17: 2 3 1 0

    index 18: 3 0 1 2
    index 19: 3 0 2 1
    index 20: 3 1 0 2
    index 21: 3 1 2 0
    index 22: 3 2 0 1
    index 23: 3 2 1 0

    --- nth Premutation:    --- get m = length and n = the index of the needed permutation (0-based)
                                let array perm of length m, numbers of length m and has numbers from 0 to m-1
                                for all numbers i form m-1 to 1
                                    let p = n / i! [index of the next number in the remaining numbers]
                                    perm[m - 1 - p] = numbers[p]
                                    erase index p from numbers
                                    n := n mod i!
                                return perm

    --- Index of permutaion:    --- let perm a the permutation of length m
                                    let idx=0
                                    for all number i from 0 to m-1
                                        idx += (n-1-i)! * prem[i]
                                        for all numbers j from i+1 to m-1
                                            if perm[j] > perm[i], prem[j] := perm[j] - 1
                                    return idx

    --- Permutation Multiplication(Mapping):    --- mapping a permutaion a = [a1, a2, .., ak] to a permutaion b = [b1, b2, .., bk]
                                                    means let a = [a[b1], a[b2], .., a[bk]] and is denoted by a*b
                                                --- It is associative, and not commutative
                                                    a*b*c = (a*b)*c = a*(b*c)
                                                    a*b != b*a
                                                --- We can use the same logic of fast power to mapping a permutation a huge number of times
                                                    a^12 = a^6 * a^6 (Order of log(power))
                                                --- Can we do better? YES! using permutation cycles

    --- Permutation Cycles: --- Cycle: Each number maps to a number, maps to a number, ..., and return to the first number
                                1 -> 3, 3 -> 4, 4 -> 2, 2 -> 1
                                ex: 2 0 1 4 3
                                    0 -> 2
                                    2 -> 1
                                    1 -> 0 (Odd Cylce)

                                    3 -> 4
                                    4 -> 3 (Even Cycle)
                            --- A cycle of length n, if we applied n time , it backs to its origin
                                n+1 time -> 1 time
                                n+2 time -> 2 times
                                n+3 time -> 3 times
                                ...
                                m times  -> m%n times
                            --- If a permutation perm has x cycles and we apply perm on itself
                                odd cycles will remains, and even cycles will be divided to 2 cycle with half of the size
                                and the elements of the original disjoint cycles will never be mixed
    --- Permutaion Order:   --- Remember: Each cycle with length ni repeats itself every ni times
                            --- So, if you have x cycles with length n1, n2, .., nx, the permutation will repeat itself
                                every lcm(n1, n2, .., nx) times
    --- Stirleng Numbers:   Again? WHO CARES???? OK, we will back to him
    --- The video will discuss some problems, please listen to them
*/

/// ========================================================================================================================================================================================================

vector<int> nth_permutation(int len, int idx)
{
    vector<int> identity(len), ans(len, 0);
    for (int i = 0; i < len; i++)
        identity[i] = i;

    for (int i = 0; i < len; i++)
    {
        int p = idx / fact[len - 1 - i];
        ans[i] = identity[p];
        identity.erase(identity.begin() + p);
        idx %= fact[len - 1 - i];
    }
    return ans;
}
int permutation_index(vector<int> perm)
{
    int n = perm.size(), ans = 0;

    for (int i = 0; i < perm.size(); i++)
    {
        ans += fact[n - 1 - i] * perm[i];
        for (int j = i + 1; j < n; j++)
            perm[j] -= (perm[j] > perm[i]);
    }
    return ans;
}
vector<int> map_permutation(vector<int> perm, vector<int> mapping)
{
    int n = perm.size();
    vector<int> ans(n);

    for (int i = 0; i < n; i++)
        ans[i] = perm[mapping[i]];
    return ans;
}
vector<vector<int>> get_permutation_cycles(vector<int> perm)
{
    int n = perm.size();
    vector<int> vis(n, 0);
    vector<vector<int>> ans;

    for (int i = 0; i < n; i++)
    {
        if (vis[i])
            continue;

        int d = perm[i];
        vector<int> temp;
        while (!vis[d])
        {
            temp.emplace_back(d);
            vis[d] = 1, d = perm[d];
        }

        ans.emplace_back(temp);
    }
    return ans;
}
int permutation_order(vector<int> perm)
{
    int ans = 1;
    vector<vector<int>> x = get_permutation_cycles(perm);

    for (auto v : x)
        ans = lcm(ans, v.size());
    return ans;
}

int32_t main()
{
    fast();
    all_prime();
    all_fact();
    all_phi();
    all_mobius();
    all_mod_inv(1e9 + 7);

    return 0;
}