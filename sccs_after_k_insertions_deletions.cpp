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

class HeavyPathDecomposition
{
private:
    vector<vector<int>> adj; // original graph adjacency list
    // vector<int> parent;            // parent in DFS tree
    vector<int> depth;    // depth in DFS tree
    vector<int> size;     // subtree size
    vector<int> heavy;    // heavy child (-1 if leaf)
    vector<bool> visited; // visited nodes during DFS

public:
    map<int, vector<vector<int>>> paths; // resulting heavy paths
    vector<vector<int>> tree_adj;        // DFS tree adjacency list (only tree edges)

    HeavyPathDecomposition(const vector<vector<int>> &graph) : adj(graph)
    {
        int n = adj.size();
        tree_adj.resize(n);
        // parent.resize(n, -1);
        depth.resize(n, 0);
        size.resize(n, 0);
        heavy.resize(n, -1);
        visited.resize(n, false);

        getHeavyPaths();
    }

    // First DFS pass to build the DFS tree and compute size/depth/heavy children
    void buildDFSTree(int u)
    {
        visited[u] = true;
        size[u] = 1;
        int max_size = 0;

        for (int v : adj[u])
        {
            if (!visited[v])
            {
                // This is a tree edge
                tree_adj[u].push_back(v);
                // parent[v] = u;
                depth[v] = depth[u] + 1;
                buildDFSTree(v);
                size[u] += size[v];
                if (size[v] > max_size)
                {
                    max_size = size[v];
                    heavy[u] = v;
                }
            }
            // Ignore back edges and cross edges
        }
    }

    // Second pass to decompose into heavy paths
    void decompose(int u, vector<int> &current_path)
    {
        current_path.push_back(u);

        // If u has a heavy child, continue the current path
        if (heavy[u] != -1)
        {
            decompose(heavy[u], current_path);
        }

        // Process light children (start new paths)
        for (int v : tree_adj[u])
        {
            if (v != heavy[u])
            { // This is a light edge
                vector<int> new_path;
                decompose(v, new_path);
                if (!new_path.empty())
                {
                    paths[depth[new_path[0]]].push_back(new_path);
                }
            }
        }
    }

    void getHeavyPaths()
    {
        if (adj.empty())
            return;

        // Build DFS tree starting from root (0)
        buildDFSTree(0);

        // Perform heavy path decomposition
        vector<int> initial_path;
        decompose(0, initial_path);
        if (!initial_path.empty())
        {
            paths[depth[initial_path[0]]].push_back(initial_path);
        }

        // return paths;
    }
};

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
map<int, vector<int>> removeFailedEdgesMap(map<int, vector<int>> &adj)
{
    map<int, vector<int>> new_adj = adj;
    for (const auto &edge : F)
    {
        int u = edge.first;
        int v = edge.second;

        // Find and remove v from u's adjacency list
        auto it = find(new_adj[u].begin(), new_adj[u].end(), v);
        if (it != new_adj[u].end())
        {
            new_adj[u].erase(it);
        }
    }
    return new_adj;
}

set<int> reachableVertices2D(int x, vector<vector<int>> &adj)
{
    set<int> B;
    B.insert(x);
    queue<int> q;
    vector<int> vis(n, 0);
    vis[x] = 1;
    q.push(x);
    while (!q.empty())
    {
        auto u = q.front();
        q.pop();
        for (auto v : adj[u])
        {
            if (!vis[v])
            {
                vis[v] = 1;
                B.insert(v);
                q.push(v);
            }
        }
    }
    return B;
}

set<int> reachableVertices(int x, map<int, vector<int>> &adj)
{
    set<int> B;

    queue<int> q;
    vector<int> vis(n, 0);
    vis[x] = 1;
    q.push(x);
    while (!q.empty())
    {
        auto u = q.front();
        q.pop();
        for (auto v : adj[u])
        {
            if (!vis[v])
            {
                vis[v] = 1;
                B.insert(v);
                q.push(v);
            }
        }
    }
    return B;
}

map<int, vector<int>> inducedSubgraph(const vector<vector<int>> &g, set<int> &A)
{
    map<int, vector<int>> sg;
    for (auto u : A)
    {
        for (auto v : g[u])
        {
            if (A.find(v) != A.end())
                sg[u].push_back(v);
        }
    }
    return sg;
}

void printMapAdj(const map<int, vector<int>> &adj)
{
    for (auto x : adj)
    {
        cout << x.first << ": ";
        for (auto v : x.second)
            cout << v << " ";
        cout << endl;
    }
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
// vertices reachable from x in set A | check by using subgraph of kftrs[x] induced by vertices in A | Lemma 4 & 10
set<int> reach(int x, set<int> &A, bool isrev)
{
    vector<vector<int>> H = isrev ? reverse_kftrs[x] : kftrs[x];

    map<int, vector<int>> g_a = inducedSubgraph(removeFailedEdges(H), A);
    printMapAdj(g_a);

    set<int> B = reachableVertices(x, g_a);
    if (A.find(x) != A.end())
        B.insert(x); // if x = path[mid] is in set A then it's reachable from itself (we r assuming this)
    return B;
}

// finds xout for vertices in A which is in path[i...j]
void binarySearch(int i, int j, set<int> &A, vector<int> &Xinfo, vector<int> &path, bool isrev)
{
    if (A.empty() || i > j)
        return;

    if (i + 1 == j)
    {
        if (isrev)
        {
            auto S = reach(path[i], A, isrev); // vertices reachable from path[i];
            for (auto x : S)
                Xinfo[x] = i;
            for (auto u : A)
            {
                if (S.find(u) == S.end())
                    Xinfo[u] = j;
            }
        }
        else
        {
            auto S = reach(path[j], A, isrev); // vertices reachable from path[i];
            for (auto x : S)
                Xinfo[x] = j;
            for (auto u : A)
            {
                if (S.find(u) == S.end())
                    Xinfo[u] = i;
            }
        }
        return;
    }
    if (i == j)
    {
        for (auto v : A)
            Xinfo[v] = i;
    }

    else
    {
        int mid = (i + j + 1) / 2; // ceil of i+j/2
        set<int> B = reach(path[mid], A, isrev);

        cout << i << " " << j << " : set B (reachable from ) " << path[mid] << " in G(A)-F\n";
        for (auto u : B)
            cout << u << " ";
        cout << endl;
        set<int> A_minus_B;
        for (auto u : A)
        {
            if (B.find(u) == B.end())
                A_minus_B.insert(u);
        }
        if (isrev)
        {
            binarySearch(mid + 1, j, A_minus_B, Xinfo, path, isrev);
            binarySearch(i, mid, B, Xinfo, path, isrev);
        }
        else
        {
            binarySearch(i, mid - 1, A_minus_B, Xinfo, path, isrev);
            binarySearch(mid, j, B, Xinfo, path, isrev);
        }
    }
}

void X_Out(vector<int> &Xout, vector<int> &X, set<int> &A)
{
    // get Xout and Xin of set A using binary search
    int t = X.size();

    //
    cout << "X[0]: " << X[0] << '\n';
    cout << "kftrs X[0] :\n";
    printAdj(kftrs[X[0]]);
    cout << "kftrsX[0]-F :\n";
    printAdj(removeFailedEdges(kftrs[X[0]]));
    cout << "Induced subgraph of kftrsX[0]-F for set A :\n";
    map<int, vector<int>> gx1 = inducedSubgraph(removeFailedEdges(kftrs[X[0]]), A);
    printMapAdj(gx1);
    ////

    cout << "X[t-1]: " << X[t - 1] << '\n';
    cout << "kftrs X[t-1] :\n";
    printAdj(kftrs[X[t - 1]]);
    cout << "kftrsX[t-1]-F :\n";
    printAdj(removeFailedEdges(kftrs[X[t - 1]]));
    cout << "Induced subgraph of kftrsX[t-1]-F for set A :\n";
    map<int, vector<int>> gxt = inducedSubgraph(removeFailedEdges(kftrs[X[t - 1]]), A);
    printMapAdj(gxt);

    // map<int, vector<int>> gx1 = inducedSubgraph(removeFailedEdges(kftrs[X[0]]), A);
    // map<int, vector<int>> gxt = inducedSubgraph(removeFailedEdges(kftrs[X[t - 1]]), A);
    // cout << "Gxt : \n";
    // printMapAdj(gxt);

    // // vertices reachable from X[t-1] in graph G[xt](A)\F | induced subgraph A of kftrs
    set<int> Vt = reachableVertices(X[t - 1], gxt); // since checking reachability from X[t] , we'll use kftrs(t)
    Vt.insert(X[t - 1]);
    cout << "Computed Vt : ";
    for (auto a : Vt)
        cout << a << " ";
    cout << endl;
    set<int> V1 = reachableVertices(X[0], gx1);
    V1.insert(X[0]);
    cout << "Computed V1 : ";
    for (auto a : V1)
        cout << a << " ";
    cout << endl;
    set<int> S; // V1\Vt;
    cout << "Computed S = V1-Vt : ";
    for (auto u : V1)
    {
        if (Vt.find(u) == Vt.end())
        {
            S.insert(u);
            cout << u << " ";
        }
    }
    cout << endl;
    for (auto v : Vt)
    {
        Xout[v] = t - 1;
    }

    binarySearch(0, t - 2, S, Xout, X, 0);
}
void X_in(vector<int> &Xin, vector<int> &X, set<int> &A)
{
    // get Xout and Xin of set A using binary search
    int t = X.size();
    // map<int, vector<int>> gx1 = inducedSubgraph(removeFailedEdges(reverse_kftrs[X[0]]), A);
    // map<int, vector<int>> gxt = inducedSubgraph(removeFailedEdges(reverse_kftrs[X[t - 1]]), A);

    cout << '\n';
    cout << "Calculation done for Xin : \n";
    cout << "X[0]: " << X[0] << '\n';
    cout << "reverse kftrs X[0] :\n";
    printAdj(reverse_kftrs[X[0]]);
    cout << "reverse kftrsX[0]-F :\n";
    printAdj(removeFailedEdges(reverse_kftrs[X[0]], true));
    cout << "Induced subgraph of reverse kftrsX[0]-F for set A :\n";
    map<int, vector<int>> gx1 = inducedSubgraph(removeFailedEdges(reverse_kftrs[X[0]], true), A);
    printMapAdj(gx1);
    ////

    cout << "X[t-1]: " << X[t - 1] << '\n';
    cout << "reverse_kftrs X[t-1] :\n";
    printAdj(reverse_kftrs[X[t - 1]]);
    cout << "reverse kftrsX[t-1]-F :\n";
    printAdj(removeFailedEdges(reverse_kftrs[X[t - 1]], true));
    cout << "Induced subgraph of reverse kftrsX[t-1]-F for set A :\n";
    map<int, vector<int>> gxt = inducedSubgraph(removeFailedEdges(reverse_kftrs[X[t - 1]], true), A);
    printMapAdj(gxt);

    // vertices reachable from X[t-1] in graph G[xt](A)\F | induced subgraph A of kftrs
    set<int> Vt = reachableVertices(X[t - 1], gxt); // since checking reachability from X[t] , we'll use kftrs(t)
    Vt.insert(X[t - 1]);
    set<int> V1 = reachableVertices(X[0], gx1);
    V1.insert(X[0]);

    cout << "Computed Vt : ";
    for (auto a : Vt)
        cout << a << " ";
    cout << endl;

    cout << "Computed V1 : ";
    for (auto a : V1)
        cout << a << " ";
    cout << endl;

    cout << "Computed S = Vt-V1 : ";

    set<int> S; // V1\Vt;
    for (auto u : Vt)
    {
        if (V1.find(u) == V1.end())
        {
            S.insert(u);
            cout << u << " ";
        }
        cout << '\n';
    }
    for (auto v : V1)
    {
        Xin[v] = 0;
    }

    binarySearch(1, t - 1, S, Xin, X, 1);
}

// SCC's of  G(A)\F intersecting path(a,b)=X |  subgraph induced by vertices in T(a) = A i.e, G(A)\F ,
vector<set<int>> get_Scc(vector<int> &X, set<int> &A)
{
    vector<int> Xout(n, -1), Xin(n, -1);
    X_Out(Xout, X, A);
    cout << "Computed Xout: \n";

    for (auto a : Xout)
        cout << a << " ";
    cout << endl;
    cout << "Computing Xin: \n";
    X_in(Xin, X, A); // these set Xin and Xout
    for (auto a : Xin)
        cout << a << " ";
    cout << endl;

    map<pair<int, int>, vector<int>> mp;
    for (int i = 0; i < n; i++)
    {
        if (Xin[i] != -1 && Xin[i] <= Xout[i])
            mp[make_pair(Xin[i], Xout[i])].push_back(i);
        // cout << i << ": " << Xin[i] << " " << Xout[i] << '\n';
    }

    cout << "Map {Xin, Xout}--> comp :\n";
    for (auto x : mp)
    {
        cout << x.first.first << " " << x.first.second << ":  ";
        for (auto u : x.second)
            cout << u << " ";
        cout << '\n';
    }

    vector<set<int>> scc;
    for (auto &x : mp)
    {
        if (!x.second.empty())
        {
            set<int> comp(x.second.begin(), x.second.end());
            scc.push_back(comp);
        }
    }
    return scc;
}

set<set<int>> solve(Graph &G)
{
    cout << "Entering solve:\n";
    set<set<int>> S; // This will be SCC's in G\F
    set<int> W;      // W is subset of V whose SCC's have been computed

    cout << "Starting hpd...\n";
    HeavyPathDecomposition hpd(G.original_adj);
    auto paths = hpd.paths; // heavy path decomposition of a dfs tree of G and paths stored in increasing order of their levels
    for (auto p : paths)
    {
        for (auto x : p.second)
        {
            for (auto u : x)
                cout << u << " ";
            cout << '\n';
        }
    }
    cout << "Hpd Complete!\n";

    for (auto &path : paths)
    {

        for (auto &X : path.second) // same depth can have multiple paths
        {
            vector<vector<int>> subpaths;
            vector<int> p;
            p.push_back(X[0]);
            for (int i = 1; i < X.size(); i++) // break a path in multiple in case it 's not intact completely
            {
                if (!inF[make_pair(X[i - 1], X[i])]) // not in F
                    p.push_back(X[i]);
                else
                {
                    subpaths.push_back(p);
                    p = {X[i]};
                }
            }
            subpaths.push_back(p);

            // compute scc intersecting p, for each path in G\F
            for (auto p : subpaths)
            {
                // for (auto u : p)
                // {
                //     cout << u << " ";
                // }
                // cout << '\n';
                cout << "\n\n";

                set<int> A = reachableVertices2D(p[0], hpd.tree_adj); // vertices lying in this subtree of dfs tree from 0
                A.insert(p[0]);                                       // p[0] is also part of subtree
                // all scc we get "now" which intersect p,  lie completely in T(a) = A , So instead of processing G*\F , we'll use kftrs(G(A))\F

                cout << "Computing SCC for path p: ";
                for (auto u : p)
                {
                    cout << u << " ";
                }
                cout << '\n';
                auto scc = get_Scc(p, A); // scc intersecting p which completely lie in T(a)
                cout << "SCC of path p done! Number of comps intersecting this path : " << scc.size() << '\n';
                // since some scc in T(a) may be part of a bigger scc as we are dealing with G*(A)\F hence
                for (auto comp : scc)
                {
                    // if  comp subset of w => comp intersect w  (disjoint sccs and earlier computed are bigger as we are dealing with smaller graph here)
                    // so if two intersect then  one is subset and has already been included so we skip that

                    bool isSubset = false;
                    cout << "{ ";
                    for (auto u : comp)
                    {
                        if (W.find(u) != W.end())
                        {
                            isSubset = true;
                            break;
                        }
                        cout << u << " ";
                    }
                    cout << " }\n";
                    if (!isSubset)
                    {
                        S.insert(comp);
                        W.insert(comp.begin(), comp.end());
                    }
                    cout << "W : ";
                    for (auto x : W)
                        cout << x << " ";
                    cout << '\n';
                }
            }
        }
    }
    return S;
}

vector<vector<int>> H;
vector<vector<int>> H_rev; // Reversed graph
vector<bool> visited;
stack<int> finish_order;

void dfs1(int u)
{
    visited[u] = true;
    for (int v : H[u])
    {
        if (!visited[v])
            dfs1(v);
    }
    finish_order.push(u);
}

void dfs2(int u, set<int> &component)
{
    visited[u] = true;
    component.insert(u);
    for (int v : H_rev[u])
    {
        if (!visited[v])
            dfs2(v, component);
    }
}

set<set<int>> kosaraju()
{
    set<set<int>> sccs;

    // Step 1: First DFS on original graph
    visited.assign(n, false);
    for (int i = 0; i < n; ++i)
    {
        if (!visited[i])
            dfs1(i);
    }

    // Step 2: Reverse the graph
    H_rev.assign(n, {});
    for (int u = 0; u < n; ++u)
    {
        for (int v : H[u])
        {
            H_rev[v].push_back(u);
        }
    }

    // Step 3: Second DFS on reversed graph in decreasing order of finish times
    visited.assign(n, false);
    while (!finish_order.empty())
    {
        int u = finish_order.top();
        finish_order.pop();
        if (!visited[u])
        {
            set<int> component;
            dfs2(u, component);
            sccs.insert(component);
        }
    }

    return sccs;
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

set<set<int>> updateSCCs(set<set<int>> &C1, set<set<int>> &C2, const unordered_set<int> &S)
{
    set<set<int>> sccs;
    for (int v : S)
    {
        set<int> a;
        // Find SCC in C2 that contains v
        for (const auto &scc : C2)
        {
            if (scc.count(v))
            {
                a = scc;
                break;
            }
        }
        if (a.empty())
            continue; // Skip if v not found in any SCC

        // Collect SCCs from C1 that are subsets of a
        vector<set<int>> toRemove;
        for (const auto &b : C1)
        {
            bool isSubset = true;
            for (int node : b)
            {
                if (!a.count(node))
                {
                    isSubset = false;
                    break;
                }
            }
            if (isSubset)
            {
                toRemove.push_back(b);
            }
        }

        // Remove those from C1
        for (const auto &b : toRemove)
        {
            C1.erase(b);
        }
        sccs.insert(a);
    }

    // Insert remaining SCCs in C1 into C2
    for (const auto &b : C1)
    {
        sccs.insert(b);
    }
    return sccs;
}

set<set<int>> fun(Graph &G)
{
    auto C1 = solve(G); // C is set of sccs in G\F.

    unordered_set<int> S;
    for (auto e : I)
    {
        S.insert(e.first);
        S.insert(e.second); // endpoints of edges in I
    }

    // Add k-FTRS_out and k-FTRS_in for vertices in S, excluding edges in X
    computeH(S);
    cout << "\n\n H : \n";
    // printAdj(H);
    H = removeFailedEdges(H);
    printAdj(H);

    auto C2 = kosaraju();
    int i = 0;
    cout << "C2\n";
    for (auto s : C2)
    {
        cout << ++i << ": { ";
        for (auto u : s)
            cout << u << " ";
        cout << "}\n";
    }

    i = 0;
    cout << "C1\n";
    for (auto s : C1)
    {
        cout << ++i << ": { ";
        for (auto u : s)
            cout << u << " ";
        cout << "}\n";
    }

    // update sccs
    auto ans = updateSCCs(C1, C2, S); // C2 will have all sccs

    return ans;
}

int main()
{
    // Take input
    // IOS;

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
        // inF[make_pair(u, v)] = 1;
    }

      auto sccs = fun(G);

    cout << "\n\nSSCs after updating edges of Graph G: \n";
    i = 0;
    for (auto s : sccs)
    {
        cout << ++i << ": { ";
        for (auto u : s)
            cout << u << " ";
        cout << "}\n";
    }
    return 0;
}