#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

int main(int argc, char **argv) {
    cin.tie(0);
    ios_base::sync_with_stdio(false);
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int t, n, m, D;
    vector<vector<int>> A, K, new_A;
    string filename;

    if (rank == 0) {
        cin >> filename;

        ifstream file(filename);

        file >> t;
        file >> n >> m;
        A.resize(n, vector<int>(m));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                file >> A[i][j];
            }
        }
        file >> D;
        K.resize(D, vector<int>(D));
        for (int i = 0; i < D; i++) {
            for (int j = 0; j < D; j++) {
                file >> K[i][j];
            }
        }

        file.close();
    }

    MPI_Bcast(&t, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        K.resize(D, vector<int>(D));
    }
    for (int i = 0; i < D; i++) {
        MPI_Bcast(K[i].data(), D, MPI_INT, 0, MPI_COMM_WORLD);
    }

    int rowsPerProc = n / size;
    int padding = (D-1) / 2;
    int startRow = rank * rowsPerProc;
    int endRow = (rank == size - 1) ? n : startRow + rowsPerProc;

    vector<vector<int>> local_A(endRow - startRow + 2 * padding, vector<int>(m, 0));
    new_A.resize(endRow - startRow + 2 * padding, vector<int>(m, 0));

    // Resize the matrix A with padding
    // A.resize(n + 2 * padding, vector<int>(m, 0));

    if (rank == 0) {
        // Distribute rows to other processes
        for (int r = 1; r < size; r++) {
            int sendStartRow = r * rowsPerProc;
            int sendEndRow = (r == size - 1) ? n : (r + 1) * rowsPerProc;
            // cout << "Rank " << r << ": " << sendStartRow << " " << sendEndRow << "\n";

            for (int i = sendStartRow; i < sendEndRow; i++) {
                MPI_Send(A[i].data(), m, MPI_INT, r, 0, MPI_COMM_WORLD);
            }
        }

        // assign rows to local_A
        for (int i = 0; i < endRow - startRow; i++) {
            local_A[i + padding] = A[i];
        }
    } else {
        // Receive rows
        for (int i = 0; i < endRow - startRow; i++) {
            MPI_Recv(local_A[i + padding].data(), m, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    // Test: print local_A with rank order
    // for (int r = 0; r < size; r++) {
    //     if (rank == r) {
    //         cout << "Rank " << rank << ":\n";
    //         for (int i = 0; i < local_A.size(); i++) {
    //             for (int j = 0; j < m; j++) {
    //                 cout << local_A[i][j] << " ";
    //             }
    //             cout << "\n";
    //         }
    //     }
    //     MPI_Barrier(MPI_COMM_WORLD);
    // }

    for (int time = 0; time < t; time++) {
        vector<int> uppersend, lowersend, upperrecv(padding * m), lowerrecv(padding * m);
        for (int i = padding; i < 2 * padding; i++) {
            for (int j = 0; j < m; j++) {
                uppersend.push_back(local_A[i][j]);
            }
        }

        // size: 4
        // [2 - 3]
        for (int i = local_A.size() - padding-1; i < local_A.size(); i++) {
            for (int j = 0; j < m; j++) {
                lowersend.push_back(local_A[i][j]);
            }
        }

        // test send size
        // cout << "Rank " << rank << ":\n";
        // cout << "m: " << m << ", padding: " << padding << "\n";
        // cout << "uppersend: " << uppersend.size() << "\n";
        // cout << "lowersend: " << lowersend.size() << "\n";

        MPI_Send(uppersend.data(), m * padding, MPI_INT, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);
        MPI_Send(lowersend.data(), m * padding, MPI_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);

        MPI_Recv(upperrecv.data(), m * padding, MPI_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(lowerrecv.data(), m * padding, MPI_INT, (rank + 1) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < padding; i++) {
            for (int j = 0; j < m; j++) {
                local_A[i][j] = upperrecv[(i)*m + j];
            }
        }
        for (int i = 0; i < padding; i++) {
            for (int j = 0; j < m; j++) {
                local_A[local_A.size() - padding + i][j] = lowerrecv[i * m + j];
            }
        }

        for (int i = padding; i < local_A.size() - padding; i++) {
            for (int j = 0; j < m; j++) {
                int sum = 0;
                for (int di = -padding; di <= padding; di++) {
                    for (int dj = -padding; dj <= padding; dj++) {
                        int x = i + di;
                        int y = (j + dj + m) % m;
                        sum += K[di + padding][dj + padding] * local_A[x][y];
                    }
                }
                new_A[i][j] = sum / (D * D) + 0.5;
            }
        }

        local_A = new_A;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // Collect results
    if (rank == 0) {
        for (int r = 1; r < size; r++) {
            int recvStartRow = r * rowsPerProc;
            int recvEndRow = (r == size - 1) ? n : (r + 1) * rowsPerProc;
            for (int i = recvStartRow; i < recvEndRow; i++) {
                MPI_Recv(A[i].data(), m, MPI_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for (int i = 0; i < endRow - startRow; i++) {
            A[i] = local_A[i + padding];
        }
    } else {
        for (int i = 0; i < endRow - startRow; i++) {
            MPI_Send(local_A[i + padding].data(), m, MPI_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                cout << A[i][j] << " ";
            }
            // cout << "\n";
        }
        // cout << "\n";
    }

    MPI_Finalize();
    return 0;
}
