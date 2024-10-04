// #pragma GCC optimize("O3,unroll-loops")
// #pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <climits>
#include <fstream>
#include <iostream>
#include <mpi.h>
#include <queue>
#include <vector>
using namespace std;
// #define endl '\n'
#define int long long

struct Edge {
    // int u; // from
    int32_t v; // to
    int32_t w; // weight
};

vector<Edge> G[50000];

priority_queue<pair<int, int32_t>, vector<pair<int, int32_t>>, greater<pair<int, int32_t>>> pq;

int32_t main(int32_t argc, char **argv) {
    // ios_base::sync_with_stdio(false);
    // cin.tie(NULL);
    MPI_Init(&argc, &argv);
    int32_t rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (rank == 0) {
        string filename;
        int n;
        cin >> filename;
        // cout << filename << endl;

        ifstream file(filename);
        file >> n;
        // cout << "n: " << n << endl;

        // G.resize(n);

        // read file
        // cout << "read file" << endl;
        Edge e;
        int32_t u;
        while (file >> u >> e.v >> e.w) {
            G[u].push_back(e);
            // cout << u << " " << e.v << " " << e.w << endl;
        }

        // file.close();
        // dijkstra
        vector<int> dist(n);
        vector<bool> vis(n);
        for (int i = 0; i < n; i++) {
            dist[i] = 1000000000;
            vis[i] = false;
        }
        dist[0] = 0;

        
        pq.push({0, 0});

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();
            if (vis[u]) {
                continue;
            }
            vis[u] = true;
            for (const auto &e : G[u]) {
                int v = e.v;
                int w = e.w;
                if (dist[v] > dist[u] + w) {
                    dist[v] = dist[u] + w;
                    pq.push({dist[v], v});
                }
            }
        }

        for (int i = 0; i < n; i++) {
            cout << dist[i] << " ";
        }
        // cout << endl;
    }
    return 0;
}