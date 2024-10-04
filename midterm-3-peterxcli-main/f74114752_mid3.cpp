#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

int mod(int a, int b) {
    return (a % b + b) % b;
}

// merge and sort
vector<int> merge_arr(const vector<int> &left_arr, const vector<int> &right_arr) {
    vector<int> new_arr;
    for (int i = 0; i < left_arr.size(); i++) new_arr.emplace_back(left_arr[i]);
    for (int i = 0; i < right_arr.size(); i++) new_arr.emplace_back(right_arr[i]);
    sort(new_arr.begin(), new_arr.end());

    return new_arr;
}

vector<int> a;
vector<int> local_a;

int32_t main(int32_t argc, char **argv) {
    int32_t rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int n;
    string filename;
    
    if (rank == 0) {
        cin >> filename;
        ifstream fin(filename);
        fin >> n;
        a.resize(n);
        for (int i = 0; i < n; i++) {
            fin >> a[i];
        }
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    int local_size = n / size;
    
    vector<int> sendcounts(size, n / size);
    vector<int> displs(size, 0);
    for (int i = 0; i < size; ++i) {
        if (i < n % size) {
            sendcounts[i]++; 
        }
        if (i > 0) {
            displs[i] = displs[i - 1] + sendcounts[i - 1];
        }
    }

    local_a.resize(sendcounts[rank]);
    MPI_Scatterv(a.data(), sendcounts.data(), displs.data(),
               MPI_INT32_T, local_a.data(), sendcounts[rank],
               MPI_INT32_T, 0, MPI_COMM_WORLD);

    sort(a.begin(), a.end());

    for (int step = 1; step < size; step *= 2) {
        if (rank % (2 * step) == 0) {
            if (rank + step < size) {
                int received_size;
                MPI_Recv(&received_size, 1, MPI_INT, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                vector<int> received_arr(received_size);
                MPI_Recv(received_arr.data(), received_size, MPI_INT32_T, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                local_a = merge_arr(local_a, received_arr);
            }
        } else {
            // Odd ranked process: send convex hull to left neighbor
            int sending_size = local_a.size();
            MPI_Send(&sending_size, 1, MPI_INT, rank - step, 0, MPI_COMM_WORLD);
            MPI_Send(local_a.data(), sending_size, MPI_INT, rank - step, 0, MPI_COMM_WORLD);
            break;
        }
    }

    if (rank == 0) {
        for (const auto& x : local_a) {
            cout << x << " ";
        }
    }

    MPI_Finalize();
    return 0;
}
