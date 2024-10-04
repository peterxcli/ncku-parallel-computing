#include <bits/stdc++.h>
#include <mpi.h>
#include <omp.h>

using namespace std;

const double alpha = 1.0; // Pheromone exponent
const double beta = 5.0;  // Distance exponent
const double rho = 0.5;   // Evaporation rate
vector<vector<int>> G;
vector<vector<double>> pheromone; // Pheromone matrix
vector<int> best_path;            // Store the best path found so far
vector<int> local_best_path;
double best_path_length = numeric_limits<double>::max();
double local_best_path_length = numeric_limits<double>::max();

// Function to initialize pheromone matrix
void initializePheromone(int n) {
    pheromone.resize(n, vector<double>(n, 1.0));
}

// Function to update pheromone trails
void updatePheromones(const vector<int> &path, double path_length) {
    double pheromone_increase = 1.0 / path_length;
    for (int i = 0; i < path.size() - 1; ++i) {
        int from = path[i];
        int to = path[i + 1];
        pheromone[from][to] = (1.0 - rho) * pheromone[from][to] + pheromone_increase;
        pheromone[to][from] = (1.0 - rho) * pheromone[to][from] + pheromone_increase;
    }
}

// Function to choose the next node based on pheromone levels
int chooseNextNode(int current, const vector<bool> &visited) {
    vector<double> probabilities(G.size(), 0.0);
    double sum = 0.0;
    for (int i = 0; i < G.size(); ++i) {
        if (!visited[i]) {
            probabilities[i] = pow(pheromone[current][i], alpha) * pow(1.0 / G[current][i], beta);
            sum += probabilities[i];
        }
    }
    double random = (double)rand() / RAND_MAX;
    double cumulative_probability = 0.0;
    int highest_node = 0;
    for (int i = 0; i < G.size(); ++i) {
        if (!visited[i]) {
            if (probabilities[i] > probabilities[highest_node]) {
                highest_node = i;
            }
            cumulative_probability += probabilities[i] / sum;
            if (random <= cumulative_probability) {
                return i;
            }
        }
    }
    return highest_node;
}

double calculatePathLength(const vector<int> &path) {
    double length = 0.0;
    for (int i = 0; i < path.size() - 1; ++i) {
        length += G[path[i]][path[i + 1]];
    }
    return length;
}

// Function to simulate ant behavior
void simulateAnts(int num_ants, int n, int iter, int rank, int size) {

    int ans_per_process = num_ants / size;
    int start = rank * ans_per_process + 1; // [start, end)
    int end = rank == size - 1 ? num_ants + 1 : start + ans_per_process;

    vector<int> ants_node(num_ants + 1);
    for (int i = 1; i <= num_ants; ++i) {
        ants_node[i] = rand() % n;
    }
    for (int t = 0; t < iter; t++) {
        vector<vector<bool>> visited(num_ants + 1, vector<bool>(n, false));
        vector<vector<int>> paths(num_ants + 1, vector<int>(n, 0));
#pragma omp parallel for num_threads(8)
        for (int i = 0; i < G.size(); ++i) {
#pragma omp parallel for num_threads(8)
            for (int k = start; k < end; k++) {
                if (i == 0) {
                    visited[k][ants_node[k]] = true;
                    paths[k][0] = ants_node[k];
                    continue;
                }
                int prev = paths[k][i - 1];
                int next = chooseNextNode(prev, visited[k]);
                paths[k][i] = next;
                visited[k][next] = true;
            }
        }

#pragma omp parallel for num_threads(8)
        for (int k = start; k < end; k++) {
            double path_length = calculatePathLength(paths[k]);
#pragma omp critical
            {
                if (path_length < local_best_path_length) {
                    local_best_path_length = path_length;
                    local_best_path = paths[k];
                }
            }
        }
#pragma omp barrier
        // Gather best_path and other information from all processes
        struct {
            double value;
            int rank;
        } in, out;

        in.value = local_best_path_length;
        in.rank = rank;

        MPI_Allreduce(&in, &out, 1, MPI_DOUBLE_INT, MPI_MINLOC, MPI_COMM_WORLD);

        if (out.value < best_path_length) {
            best_path_length = out.value;
        }
        if (out.rank == rank) {
            best_path = local_best_path;
        }
        MPI_Bcast(&best_path_length, 1, MPI_DOUBLE, out.rank, MPI_COMM_WORLD);

        MPI_Barrier(MPI_COMM_WORLD);
        updatePheromones(best_path, best_path_length);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n, m, t;
    if (rank == 0) {
        cin >> n >> m >> t;
        G.resize(n, vector<int>(n, 0));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cin >> G[i][j];
            }
        }
    }

    // Broadcast problem data to all MPI processes
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    for (int i = 0; i < n; i++) {
        MPI_Bcast(&G[i][0], n, MPI_INT, 0, MPI_COMM_WORLD);
    }

    // Initialize pheromone matrix
    initializePheromone(n);

    // Call the ant simulation function
    simulateAnts(m, n, t, rank, size);

    // Print the best path found
    if (rank == 0) {
        cout << "Best path: ";
        for (int i = 0; i < best_path.size(); ++i) {
            cout << best_path[i] << " ";
        }
        cout << endl;
        cout << "Best path length: " << best_path_length << endl;
    }

    MPI_Finalize();
    return 0;
}
