// Worked well for all testcases mentioned in OneNote

#include <bits/stdc++.h>
using namespace std;

#define IOS                  \
    ios::sync_with_stdio(0); \
    cin.tie(0);              \
    cout.tie(0);

struct Edge
{
    int v;         // Target vertex
    int rev_index; // Index of reverse edge
    bool active;   // Represents residual capacity (true = available)

    Edge(int v, int rev_index, bool is_active)
        : v(v), rev_index(rev_index), active(is_active) {}
};

class Graph
{
public:
    int n;                                    // Number of vertices
    vector<vector<int>> original_adj;         // Original graph (never modified)
    vector<vector<Edge>> residual_adj, r_adj; // Residual graph (updated during flow)
    vector<vector<int>> incoming;

    Graph(int V) : n(V), original_adj(V), residual_adj(V), incoming(V) {}

    // Add edge to both original and residual graphs   | Time: O(1) amortized (vector push_back)
    void addEdge(int u, int v)
    {
        original_adj[u].push_back(v);

        residual_adj[u].push_back(Edge(v, residual_adj[v].size(), 1));
        residual_adj[v].push_back(Edge(u, residual_adj[u].size() - 1, 0));

        incoming[v].push_back(u);
    }

    // multisource bfs to find if there exists an augmenting path in residual graph from any of the source in S to sink t
    bool bfs(set<int> &S, int t, vector<int> &parent, vector<Edge *> &parent_edge)
    {
        parent.assign(n, -1);
        parent_edge.assign(n, nullptr); // Initialize edge pointers
        queue<int> q;
        for (auto s : S)
        {
            q.push(s);
            parent[s] = -2; // we never come to s
        }

        while (!q.empty())
        {
            auto u = q.front();
            q.pop();

            for (Edge &e : r_adj[u])
            {
                if (parent[e.v] == -1 && e.active)
                {
                    parent[e.v] = u;
                    parent_edge[e.v] = &e;
                    if (e.v == t)
                        return true;
                    q.push(e.v);
                }
            }
        }
        return false;
    }

    // find maxflow from S to t and update residual graph | Edmond -karp BFS
    int get_maxflow(set<int> &S, int t)
    {
        r_adj = residual_adj;
        int max_flow = 0;
        vector<int> parent;
        vector<Edge *> parent_edge;

        while (bfs(S, t, parent, parent_edge))
        { // if path exists then get that path using parent array to modify the residual graph
            int v = t;
            while (parent[v] != -2)
            {
                Edge *e = parent_edge[v];
                e->active = 0;
                r_adj[e->v][e->rev_index].active = 1;
                v = parent[v];
                // cout << v << " ";
            }
            // cout << endl;
            max_flow++;
        }
        return max_flow;
    }

    set<int> findReachableToT(int t)
    {
        set<int> B;
        queue<int> q;
        q.push(t);
        B.insert(t);

        while (!q.empty())
        {
            int u = q.front();
            q.pop();

            for (Edge &e : r_adj[u])
            {
                if (B.find(e.v) == B.end() && !e.active)
                {
                    B.insert(e.v);
                    q.push(e.v);
                }
            }
        }
        return B;
    }

    set<int> getFlowEdgesToT(int t)
    {
        set<int> flow_edges;
        for (auto u : incoming[t])
        {
            for (Edge &e : r_adj[u])
            {
                if (e.v == t && !e.active)
                {
                    flow_edges.insert(u);
                    break;
                }
            }
        }
        return flow_edges;
    }

    Graph kftrs_t(int s, int t, int k)
    {
        set<int> S = {s};
        for (int i = 1; i <= k; i++)
        {
            get_maxflow(S, t); // this will give us residual graph using which we will find our partition A , B
            auto B = findReachableToT(t);
            set<int> A;
            for (int u = 0; u < n; u++)
            {
                if (B.find(u) == B.end())
                    A.insert(u);
            }

            // S = A union (vertices in B which have incoming edges from any vertex of A)
            S = A;
            for (auto u : A)
            {
                for (int v : original_adj[u])
                {
                    if (v != t && B.find(v) != B.end())
                        S.insert(v);
                }
            }
        }

        get_maxflow(S, t);
        Graph G_star(n);
        for (int u = 0; u < n; u++)
        {
            for (auto v : original_adj[u])
            {
                if (v != t)
                    G_star.addEdge(u, v);
            }
        }

        auto flow_edges = getFlowEdgesToT(t);
        for (int u : flow_edges)
        {
            G_star.addEdge(u, t);
        }

        return G_star;
    }

    void printRGraph()
    {
        for (int i = 0; i < n; i++)
        {
            cout << i << "---> ";
            for (auto &e : r_adj[i])
            {
                cout << e.v << ": " << e.active << " , ";
            }
            cout << endl;
        }
    }
    void printGraph()
    {
        for (int i = 0; i < original_adj.size(); i++)
        {
            cout << i << "----> ";
            for (auto x : original_adj[i])
                cout << x << " ";
            cout << endl;
        }
    }
};

map<int, vector<vector<int>>> kftrs, reverse_kftrs;
vector<pair<int, int>> F, I;
map<pair<int, int>, int> inF; // if edge (u,v) is in F
int n, m, k;

void preprocessing(Graph &G, Graph &Gr)
{
    // transform part is still left

    for (int s = 0; s < G.n; s++) // for each vertex in G compute kftrs with s as source
    {
        Graph current_G = G;
        for (int i = 0; i < G.n; i++) // Locality Lemma
        {
            current_G = current_G.kftrs_t(s, i, k);
        }
        kftrs[s] = current_G.original_adj;
    }

    for (int s = 0; s < G.n; s++) // for each vertex in reverse G, compute kftrs with s as source
    {
        Graph current_G = Gr;
        for (int i = 0; i < n; i++)
        {
            current_G = current_G.kftrs_t(s, i, k);
        }
        reverse_kftrs[s] = current_G.original_adj;
    }
}

// remove F from a graph
vector<vector<int>> removeFailedEdges(vector<vector<int>> &adj, bool rev = false)
{
    vector<vector<int>> new_adj = adj;
    for (const auto &edge : F)
    {
        int u = edge.first;
        int v = edge.second;

        if (adj[u].empty())
            continue;

        if (rev)
            swap(u, v);
        // Find and remove v from u's adjacency list
        auto it = find(new_adj[u].begin(), new_adj[u].end(), v);
        if (it != new_adj[u].end())
        {
            new_adj[u].erase(it);
        }
    }
    return new_adj;
}

void printAdj(const vector<vector<int>> &adj)
{
    for (int i = 0; i < adj.size(); i++)
    {
        cout << i << ": ";
        for (auto x : adj[i])
            cout << x << " ";
        cout << endl;
    }
}

vector<vector<int>> H, H_minus_F;

bool isReachable(int x, int y, const vector<vector<int>> &adj)
{
    vector<bool> visited(n, false);
    queue<int> q;
    q.push(x);
    visited[x] = true;
    bool x_to_y = false;
    while (!q.empty() && !x_to_y)
    {
        int u = q.front();
        q.pop();
        for (int v : adj[u])
        {
            if (v == y)
            {
                x_to_y = true;
                break;
            }
            if (!visited[v])
            {
                visited[v] = true;
                q.push(v);
            }
        }
    }
    return x_to_y;
}

void computeH(unordered_set<int> &S)
{
    H.resize(n);
    map<pair<int, int>, int> fre; // fre of edges in ftrs

    for (int v : S)
    {
        map<pair<int, int>, int> temp;
        for (int i = 0; i < n; i++)
        {
            for (auto j : kftrs[v][i])
                temp[make_pair(i, j)]++; // i --> j edge
        }

        for (int i = 0; i < n; i++)
        {
            for (auto j : kftrs[v][i])
                fre[make_pair(i, j)] = max(temp[make_pair(i, j)], fre[make_pair(i, j)]); // i --> j edge
        }

        temp.clear();
        for (int i = 0; i < n; i++)
        {
            for (auto j : reverse_kftrs[v][i])
                temp[make_pair(j, i)]++; // i --> j edge
        }

        for (int i = 0; i < n; i++)
        {
            for (auto j : reverse_kftrs[v][i])
                fre[make_pair(j, i)] = max(temp[make_pair(j, i)], fre[make_pair(j, i)]); // j--->i edge | purpose of reverse kftrs is to get vertices which can reach i
        }
    }

    for (auto e : fre)
    {
        int f = e.second;
        while (f--)
        {
            H[e.first.first].push_back(e.first.second);
        }
    }

    // insert edges in set Y
    for (auto e : I)
    {
        H[e.first].push_back(e.second);
    }
}

// Helper function to compute reachable nodes
unordered_set<int> computeReachable(int start, const vector<vector<int>> &adj)
{
    unordered_set<int> visited;
    queue<int> q;
    q.push(start);
    visited.insert(start);
    while (!q.empty())
    {
        int u = q.front();
        q.pop();
        for (int v : adj[u])
        {
            if (!visited.count(v))
            {
                visited.insert(v);
                q.push(v);
            }
        }
    }
    return visited;
}
// Transpose a graph (reverse all edges)
vector<vector<int>> transposeGraph(const vector<vector<int>> &adj)
{
    vector<vector<int>> trans(adj.size());
    for (int u = 0; u < adj.size(); ++u)
    {
        for (int v : adj[u])
        {
            trans[v].push_back(u);
        }
    }
    return trans;
}

// Optimized query function
bool queryOptimized(int x, int y)
{
    unordered_set<int> S;
    for (auto &e : I)
    {
        S.insert(e.first);
        S.insert(e.second);
    }

    // Case 1: At least one endpoint is in S
    if (S.count(x) || S.count(y))
    {
        bool x_to_y = isReachable(x, y, H_minus_F);
        bool y_to_x = isReachable(y, x, H_minus_F);
        return x_to_y && y_to_x;
    }

    // Case 2: Neither x nor y is in S

    // Check if they're strongly connected in G\F using k-FTRS
    auto kFTRS_x = removeFailedEdges(kftrs[x]);
    auto kFTRS_y = removeFailedEdges(kftrs[y]);
    bool x_to_y_G = isReachable(x, y, kFTRS_x);
    bool y_to_x_G = isReachable(y, x, kFTRS_y);

    if (x_to_y_G && y_to_x_G)
        return true;

    // Precompute reachability sets
    auto Reach_x = computeReachable(x, kFTRS_x);
    auto Reach_y = computeReachable(y, kFTRS_y);

    // Build reversed H\F for inverse reachability
    auto reversed_H = transposeGraph(H_minus_F);
    auto Reach_y_inv = computeReachable(y, reversed_H); // Nodes that can reach y in H\F
    auto Reach_x_inv = computeReachable(x, reversed_H); // Nodes that can reach x in H\F

    // Check x->y via inserted edges
    bool x_to_y_via_S = false, y_to_x_via_S = false;
    if (!x_to_y_G)
    {
        for (auto &e : I)
        {
            int u = e.first, v = e.second;
            if (Reach_x.count(u) && Reach_y_inv.count(v))
            {
                x_to_y_via_S = true;
                break;
            }
        }
    }

    // Check y->x via inserted edges
    if (!y_to_x_G)
    {
        for (auto &e : I) // iterate on e to be first edge of Y in path P
        // i.e., if there exist a path in G+U from u to v which passes through any edge of I , then P1 connects
        {
            int u = e.first, v = e.second;
            if (Reach_y.count(u) && Reach_x_inv.count(v))
            {
                y_to_x_via_S = true;
                break;
            }
        }
    }

    return (x_to_y_G || x_to_y_via_S) && (y_to_x_G || y_to_x_via_S);
}

void solve()
{
    // auto C1 = solve(G); // C is set of sccs in G\F.

    unordered_set<int> S;
    for (auto e : I)
    {
        S.insert(e.second); // endpoints of edges in I
    }

    // Construct H (k-FTRS of S + inserted edges)
    // Add k-FTRS_out and k-FTRS_in for vertices in S, excluding edges in X
    computeH(S);
    cout << "\n\n H-F : \n";
    // printAdj(H);
    H_minus_F = removeFailedEdges(H);
    printAdj(H_minus_F);

    cout << "Enter number of queries for same updates : \n";
    int q;
    cin >> q;
    while (q--)
    {
        int x, y;
        cin >> x >> y;
        if (queryOptimized(x, y))
        {
            cout << x << " and " << y << " are strongly connected after updates.\n";
        }
        else
        {
            cout << x << " and " << y << " are NOT strongly connected after updates.\n";
        }
    }
}

int main()
{

    cin >> n >> m >> k;
    Graph G(n), Gr(n);
    cout << "Input the edges of Graph\n";
    for (int i = 0; i < m; i++)
    {
        int u, v;
        cin >> u >> v;
        G.addEdge(u, v);
        Gr.addEdge(v, u);
    }

    preprocessing(G, Gr); // computes kftrs and reverse graph kftrs | space : O(2^k.n^2) | time : O(2^k.n^2); ,
    // if transformed graph then number of vertices O(m) :  be O(2^k.n.m) | and once we transform back space: O(2^k.n^2)

    cout << "Input number of failed edges and then edges:\n";
    int f;
    cin >> f;

    while (f--)
    {
        int u, v;
        cin >> u >> v;
        F.push_back({u, v});
        inF[make_pair(u, v)] = 1;
    }

    cout << "Input number of inserted edges and then edges:\n";
    int i;
    cin >> i;

    while (i--)
    {
        int u, v;
        cin >> u >> v;
        I.push_back({u, v});
    }

    solve();

    return 0;
}