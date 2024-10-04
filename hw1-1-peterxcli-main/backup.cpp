#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;
#define int long long
#define max(a, b) \
    ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })

int32_t main(int32_t argc, char **argv) {
    int32_t rank, size;
    MPI_Init(&argc, &argv);               // Initialize the MPI execution environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the calling process in the communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the size of the group associated with a communicator
    // cout << "rank: " << rank << " size: " << size << endl;
    // cout << argv[1] << endl;
    string filename;
    ifstream fin;
    int n, m;
    if (rank == 0) {
        cin >> filename;
        fin.open(filename);
        fin >> n >> m;
    }
    MPI_Bcast(&n, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
    // cout << "rank: " << rank << " n: " << n << " m: " << m << endl;

    vector<vector<int>> test(m);
    if (rank == 0) {
        for (int i = 0; i < m; i++) {
            int test_num, _;
            fin >> test_num >> _;
            test[i].resize(test_num);
            for (int j = 0; j < test_num; j++) {
                fin >> test[i][j];
                test[i][j]--;
            }
        }
    }

    for (int i = 0; i < m; i++) {
        int test_size = test[i].size();
        MPI_Bcast(&test_size, 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
        if (rank != 0) {
            test[i].resize(test_size);
        }
        for (int j = 0; j < test_size; j++) {
            MPI_Bcast(&test[i][j], 1, MPI_INT64_T, 0, MPI_COMM_WORLD);
        }
        // MPI_Bcast(test[i].data(), test_size, MPI_INT64_T, 0, MPI_COMM_WORLD);
        for (int j = 0; j < test_size; j++) {
            cout << test[i][j] << endl;
        }
    }

    int chunk = (1 << m) / size;
    int start = rank * chunk;
    int end = (rank != size - 1 ? (rank + 1) * chunk : (1 << m));
    int local_ans = 0;
    for (int i = start; i < end; i++) {
        vector<bool> program(n);
        for (int j = 0; j < n; j++) {
            if ((1 << j) & i) {
                for (int k = 0; k < test[j].size(); k++) {
                    program[test[j][k]] = 1;
                }
            }
        }
        bool flag = 1;
        for (int j = 0; j < n; j++) {
            if (program[j] == 0) {
                flag = 0;
                break;
            }
        }
        if (flag) {
            local_ans++;
        }
    }
    int global_ans;
    MPI_Reduce(&local_ans, &global_ans, 1, MPI_INT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        cout << global_ans << endl;
    }

    MPI_Finalize(); // Terminates MPI execution environment
    return 0;
}
