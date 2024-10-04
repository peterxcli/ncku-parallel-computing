#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

struct Edge {
    int u; // from
    int v; // to
    int w; // weight
};

struct Node {
    int vertex;
    int distance;

    bool operator>(const Node &other) const {
        return distance > other.distance;
    }
    bool operator<(const Node &other) const {
        return distance < other.distance;
    }
};

vector<int> dijkstra(const vector<vector<Edge>> &local_G, int source, int start, int end, int n, int rank) {
    vector<int> dist(n, INT_MAX);
    map<int, map<int, int>> weight; // weight[u][v] = w
    vector<bool> visited(n, false);

    // Priority queue to find the local minimum
    priority_queue<Node, vector<Node>, greater<Node>> pq;

    set<int> sptSet;
    for (int i = start; i < end; i++) {
        for (const auto &e : local_G[i]) {
            weight[e.u][e.v] = e.w;
        }
    }

    // init dist
    for (int i = start; i < end; i++) {
        if (weight.find(source) != weight.end() && weight[source].find(i) != weight[source].end()) {
            dist[i] = weight[source][i];
        }
        sptSet.insert(i);
        pq.push({vertex : i, distance : dist[i]});
    }

    dist[source] = 0;
    // pq.push({vertex : source, distance : 0});

    // // Debug print 1: Initial distances
    // cout << "Initial distances:" << endl;
    // for (int i = start; i < end; i++) {
    //     cout << "dist[" << i << "]: " << dist[i] << " ";
    // }
    // cout << endl;

    // if (weight.find(source) != weight.end())
    //     pq.push({vertex : source, distance : 0});

    // sptSet.erase(source);
    // weight.erase(source);

    for (int i = 0; i < n - 1; i++) {
        // cout << "start: " << start << ", end: " << end << ", sptSet: " << sptSet.size() << endl;
        int min_dist = INT_MAX, min_index = -1;

        // 1. Each process finds its local minimum distance using the priority queue
        // while (!pq.empty()) {
        //     if (sptSet.find(pq.top().vertex) == sptSet.end()) {
        //         min_dist = pq.top().distance;
        //         min_index = pq.top().vertex;
        //         pq.pop();
        //     //     break;
        //     }
        //     else {
        //         // min_dist = pq.top().distance;
        //         // min_index = pq.top().vertex;
        //         // pq.pop();
        //         break;
        //     }
        // }
        if (!pq.empty()) {
            do {
                min_dist = pq.top().distance;
                min_index = pq.top().vertex;
                pq.pop();
            } while (visited[min_index] && !pq.empty());
        }

        if (sptSet.empty()) {
            min_dist = INT_MAX;
            min_index = -1;
        }

        // // Debug print 2: Local minimum
        // cout << "Round: " << i << " Rank: " << rank << " Local minimum: vertex=" << min_index << ", distance=" << min_dist << endl;

        // 2. Reduce to get global minimum
        struct {
            int value;
            int index;
        } local_min, global_min;

        local_min.value = min_dist;
        local_min.index = min_index;

        // cout << "start: " << start << ", end: " << end << ", local_min: " << local_min.value << ", " << local_min.index << endl;
        MPI_Allreduce(&local_min, &global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        // cout << "start: " << start << ", end: " << end << ", global_min: " << global_min.value << ", " << global_min.index << endl;
        // if (global_min.value == INT_MAX) break;

        // if (rank == 0)
        // cout << "Round: " << i << " Rank: " << rank << " Global minimum: vertex=" << global_min.index << ", distance=" << global_min.value << endl;

        int u = global_min.index;
        if (u >= start && u < end) {
            // dist[u] = min(dist[u], global_min.value);
            sptSet.erase(u);
            visited[u] = true;
            // pq.push({vertex : u, distance : global_min.value});
            // pq.pop();
            // sptSet.erase(u);
        }
        else {
            pq.push({vertex : min_index, distance : min_dist});
        }

        // pq.push({vertex : u, distance : global_min.value});

        // 3. Each process updates its distance values
        for (auto v : sptSet) {
            if (weight[u].find(v) == weight[u].end()) continue;
            int w = weight[u][v];
            if (global_min.value != INT_MAX && global_min.value + w < dist[v]) {
                dist[v] = min(dist[v], global_min.value + w);
                // cout << "start: " << start << ", end: " << end << " found new smaller dist" << ", dist[" << v << "]: " << dist[v] << endl;
                // cout << "Round: " << i << " Found new smaller dist" << ", dist[" << v << "]: " << dist[v] << endl;
                pq.push({vertex : v, distance : dist[v]});
            }
        }

        // // Debug print 3: Updated distances
        // cout << "Round: " << i << " Rank: " << rank << " Updated distances:" << endl;
        // for (int i = start; i < end; i++) {
        //     cout << "dist[" << i << "]: " << dist[i] << " ";
        // }
        // cout << endl;

        // MPI_Barrier(MPI_COMM_WORLD);
    }

    return dist;
}

int main(int argc, char **argv) {
    cin.tie(0);
    ios_base::sync_with_stdio(false);
    MPI_Init(&argc, &argv);
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Datatype MPI_EDGE;
    MPI_Type_contiguous(3, MPI_INT, &MPI_EDGE);
    MPI_Type_commit(&MPI_EDGE);

    string filename;
    vector<vector<Edge>> G;
    int n;

    if (rank == 0) {
        cin >> filename;

        ifstream file(filename);
        file >> n;

        G.resize(n);
        Edge e;
        while (file >> e.u >> e.v >> e.w) {
            // G[e.u].push_back(e);
            G[e.v].push_back(e);
        }

        file.close();
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    vector<vector<Edge>> local_G;
    int basic_chunk_size = n / size;
    int remainder = n % size;

    int start = rank * basic_chunk_size + (rank < remainder ? rank : remainder);
    int end = start + basic_chunk_size + (rank < remainder ? 1 : 0);

    if (rank >= n)
        return 0;

    // test
    local_G.resize(n);

    // send partitioned graph to each process
    if (rank == 0) {
        for (int i = 1; i < size; i++) {
            int s = i * basic_chunk_size + (i < remainder ? i : remainder);
            int e = s + basic_chunk_size + (i < remainder ? 1 : 0);
            for (int j = s; j < e; j++) {
                int size = G[j].size();
                MPI_Send(&size, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
                MPI_Send(&G[j][0], G[j].size(), MPI_EDGE, i, 0, MPI_COMM_WORLD);
            }
        }
        // fill local_g for rank 0
        for (int i = start; i < end; i++) {
            local_G[i] = G[i];
        }
    } else {
        for (int i = start; i < end; i++) {
            int size;
            MPI_Recv(&size, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            local_G[i].resize(size);
            MPI_Recv(local_G[i].data(), size, MPI_EDGE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // /// pretty print the receive graph
    // for (int r = 0; r < size; r++) {
    //     // test
    //     if (rank == r) {
    //         cout << "rank: " << rank << " local_G: " << endl;
    //         for (int i = start; i < end; i++) {
    //             for (int j = 0; j < local_G[i].size(); j++) {
    //                 cout << local_G[i][j].u << " " << local_G[i][j].v << " " << local_G[i][j].w << endl;
    //             }
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    int source = 0;

    vector<int> local_dists = dijkstra(local_G, source, start, end, n, rank);
    if (rank == 0) {
        vector<int> global_dists(n);
        MPI_Gather(local_dists.data() + start, end - start, MPI_INT, global_dists.data(), end - start, MPI_INT, 0, MPI_COMM_WORLD);
        for (int i = 0; i < n; i++) {
            // cout << "Distance from " << source << " to " << i << " is " << global_dists[i] << endl;
            cout << global_dists[i] << " ";
        }
    } else {
        MPI_Gather(local_dists.data() + start, end - start, MPI_INT, nullptr, 0, MPI_INT, 0, MPI_COMM_WORLD);
    }

    MPI_Type_free(&MPI_EDGE);
    MPI_Finalize();
    return 0;
}
