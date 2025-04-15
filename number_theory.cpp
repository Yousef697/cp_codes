#include <bits/stdc++.h>

using namespace std;

const int N = 1e6 + 20;
int spf[N];
vector<int> primes;
void pre()
{
    for (int i = 1; i < N; i++)
        spf[i] = i;
    for (int i = 2; i * i < N; i++)
    {
        if (spf[i] != i)
            continue;
        for (int j = i * i; j < N; j += i)
            spf[j] = min(spf[j], i);
    }
    for (int i = 2; i < N; i++)
        if (spf[i] == i)
            primes.push_back(i);
}

int32_t main()
{
    ios::sync_with_stdio(false);
    cin.tie(nullptr), cout.tie(nullptr);

    return 0;
}