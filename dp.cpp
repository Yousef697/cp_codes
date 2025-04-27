#include <bits/stdc++.h>
#define int long long
#define all(x) x.begin(), x.end()
#define rall(x) x.rbegin(), x.rend()
#define sz(x) (int)(x).size()

using namespace std;
const double pi = acos(-1);
const int mod = 1e9 + 7;
const int N = 2e5 + 5, M = 1e2 + 2, K = 5e2 + 5;
const int MAX_WEIGHT = 10005;
const int MAX_1e2 = 105;
const int MAX_1e3 = 1005;
const int MAX_5e3 = 5005;
const int MAX_1e4 = 10005;
const int MAX_1e5 = 100005;
const int MAX_1e6 = 1000005;

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
    --- Chapter: Knapsack & Building .... (Done)
    --- Chapter: LIS & Building ......... (Done)
    --- Chapter: LCS & Building ......... (Done)
    --- Chapter: Coin Exchange .......... (Not Done)
    --- Chapter: Ranges dp .............. (Not Done)
    --- Chapter: Bitmask dp ............. (Not Done)
    --- Chapter: Counting dp ............ (Not Done)
    --- Chapter: Sums ................... (Not Done)
    --- Chapter: Strings ................ (Not Done)
*/

/// ========================================================================================================================================================================================================

/// === Chapter: Knapsack
/*
    --- The main idea is: Take or leave

    --- for the recursive code: the function depends on the index and how much weight we have [ knapsack(i, rem) ]
    --- Recurrence:     --- let: a = the answer if we leave the item = knapsack(i+1, rem)
                        --- let: b = the answer if we can take the item + value of the current item = knapsack(i+1, rem - weights[i]) + values[i]
                        --- return max(a, b)

    --- for the iterative code: we loop on all items and on all remainders
    --- Recurrences:    --- dp[i][j] = we are at the i-th item with remainder j
                        --- dp[i+1][j] = leave the current item and get the answer if we go to the next item with the same remainder
                        --- dp[i+1][j-weights[i]] = take the current item and get the answer if we go to the next item with remainder j-wights[i]
*/
int dp_knapsack[MAX_1e2][MAX_WEIGHT];
vector<int> weights(MAX_1e2), values(MAX_1e2);
vector<pair<int, int>> ans_knapsack;

/// Recursive Code, Time Comlexity: O(N * W), Memory Complexity: O(N * W)
int knapsack(int i, int rem)
{
    if (i == weights.size())
        return 0;

    int &ret = dp_knapsack[i][rem];
    if (ret != -1)
        return ret;

    ret = knapsack(i + 1, rem);
    if (rem >= weights[i])
        ret = max(ret, values[i] + knapsack(i + 1, rem - weights[i]));

    return ret;
}

void build_knapsack(int i, int rem)
{
    if (i == weights.size())
        return;

    int ret = dp_knapsack[i][rem];

    if (ret == knapsack(i + 1, rem))
        build_knapsack(i + 1, rem);

    else if (rem >= weights[i] && ret == values[i] + knapsack(i + 1, rem - weights[i]))
    {
        ans_knapsack.emplace_back(weights[i], values[i]);
        build_knapsack(i + 1, rem - weights[i]);
    }
}

/// Iterative Code, Time Comlexity: O(N * W), Memory Complexity: O(N * W)
int knapsack(int rem)
{
    for (int i = weights.size() - 1; i >= 0; i--)
        for (int j = rem; j >= 0; j--)
        {
            dp_knapsack[i][j] = dp_knapsack[i + 1][j];
            if (j >= weights[i])
                dp_knapsack[i][j] = max(dp_knapsack[i][j], dp_knapsack[i + 1][j - weights[i]] + values[i]);
        }
    return dp_knapsack[0][rem];
}

/// ========================================================================================================================================================================================================

/// === Chapter: LIS
/*
    --- The main idea: Take or leave

    --- for the recursive code: the function depends on the current index and the previous index
    --- Recurrenes:     --- let: a = the answer if we go to the next index with the same previous index [ lis(i+1, prev) ]
                        --- let: b = 1 + the answer if we can take the current index and the previous index will be the current index [ lis(i+1, i) ]
                        --- return max(a, b)

    --- for the iterative code: for each index i we look for all previous indecies
        and take the index j that dp[j] is maximized and let dp[i] = dp[j]+1
    --- Recurrences:    --- dp[i] = the answer if we let the element i the last element of the answer

    --- for the second iterative code:
        --- we have a vector "seq" that will holds the numbers
        --- for each index i, get the index minimum number in "seq" the is greater than v[i]
        --- if the index is greater than the "seq" size, add v[i] to the seq
        --- if that number is greater than v[i], we modify it to v[i] because the smaller one will give us more chioces
        --- the answer of v[i] is its index in the vector "seq"
*/
vector<int> v(MAX_5e3);

/// Recursive Code, Time Comlexity: O(N^2), Memory Complexity: O(N^2)
int dp_lis[MAX_5e3][MAX_5e3];
int lis(int i, int prev)
{
    if (i == v.size())
        return 0;

    int &ret = dp_lis[i][prev];
    if (ret != -1)
        return ret;

    ret = lis(i + 1, prev);
    if (v[i] > v[prev])
        ret = max(ret, lis(i + 1, i) + 1);

    return ret;
}

/// Iterative Code, Time Comlexity: O(N^2), Memory Complexity: O(N)
int dp_lis_[MAX_5e3];
int lis()
{
    memset(dp_lis_, 1, sizeof dp_lis_);
    for (int i = 0; i < v.size(); i++)
        for (int j = 0; j < i; j++)
            if (v[i] > v[j])
                dp_lis_[i] = max(dp_lis_[i], dp_lis_[j] + 1);
    int ans = 0;
    for (int i = 0; i < v.size(); i++)
        ans = max(ans, dp_lis_[i]);
    return ans;
}

/// Iterative Code, Time Comlexity: O(N*LOG(N)), Memory Complexity: O(N)
int dp_lis__[MAX_1e5];
pair<int, vector<int>> lis()
{
    vector<int> seq;

    for (int i = 0; i < v.size(); i++)
    {
        int ind = lower_bound(all(seq), v[i]) - seq.begin();

        if (ind == seq.size())
            seq.emplace_back(v[i]);

        if (v[i] < seq[ind])
            seq[ind] = v[i];
        dp_lis__[i] = ind + 1;
    }
    int ans = 0;
    for (int i = 0; i < v.size(); i++)
        ans = max(ans, dp_lis__[i]);

    int x = ans, prev = 1e9;
    vector<int> ans_lis;
    for (int i = v.size(); i >= 0; i--)
    {
        if (dp_lis__[i] == x && v[i] < prev)
        {
            x--;
            ans_lis.emplace_back(v[i]);
            prev = v[i];
        }
    }
    return {ans, ans_lis};
}

/// ========================================================================================================================================================================================================

/// === Chapter: LCS
/*
    --- The main idea: pointer on each string, i on s, j on t, try to move i or j to match

    --- for the recursive code: the function depends on the index i of the string s and the index j
        of the string t
    --- Recurrences:    --- let a = the answer if we move i
                            let b = the answer if we move j
                            return max(a, b)

    --- for the iterative code:
        --- we loop on each character in s and t
        --- if s[i] = t[j], there is a match, get the answer of dp_lcs[i-1][j-1] and add 1 (for the matched characters)
        --- else if the answer of the previous i with the current j is greater than the answer of the current i with
            the previous j, we let dp_lcs[i][j] = dp_lcs[i-1][j]
        --- else we let dp_lcs[i][j] = dp_lcs[i][j-1]

*/
string s, t, ans_lcs = "";
int dp_lcs[MAX_1e3][MAX_1e3];

/// Recursive Code, Time Comlexity: O(N^2), Memory Complexity: O(N^2)
int lcs(int i, int j)
{
    if (i == s.size() || j == t.size())
        return 0;

    int &ret = dp_lcs[i][j];

    if (~ret)
        return ret;

    if (s[i] == t[j])
        return 1 + lcs(i + 1, j + 1);

    int a = lcs(i + 1, j);
    int b = lcs(i, j + 1);

    return ret = max(a, b);
}

void build_lcs(int i, int j)
{
    if (i == s.size() || j == t.size())
        return;

    int ret = dp_lcs[i][j];

    if (s[i] == t[j])
    {
        ans_lcs.push_back(s[i]);
        build_lcs(i + 1, j + 1);
        return;
    }

    int a = lcs(i + 1, j);
    int b = lcs(i, j + 1);

    if (a >= b)
        build_lcs(i + 1, j);
    else
        build_lcs(i, j + 1);
}

/// Iterative Code, Time Comlexity: O(N^2), Memory Complexity: O(N^2)
int lcs()
{
    int n = s.size(), m = t.size();
    for (int i = 0; i <= n; i++)
        for (int j = 0; j <= m; j++)
        {
            if (i == 0 || j == 0)
                dp_lcs[i][j] = 0;
            else if (s[i - 1] == t[i - 1])
                dp_lcs[i][j] = 1 + dp_lcs[i - 1][j - 1];
            else if (dp_lcs[i - 1][j] >= dp_lcs[i][j - 1])
                dp_lcs[i][j] = dp_lcs[i - 1][j];
            else
                dp_lcs[i][j] = dp_lcs[i][j - 1];
        }
    return dp_lcs[n][m];
}

/// ========================================================================================================================================================================================================

int32_t main()
{
    fast();

    return 0;
}