#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;

#define int long long
int32_t main(int32_t argc, char **argv) {
    cin.tie(0);
    ios_base::sync_with_stdio(false);
    int32_t rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    int t, n, m, D1, D2;
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
        file >> D1 >> D2;
        K.resize(D1, vector<int>(D2));
        for (int i = 0; i < D1; i++) {
            for (int j = 0; j < D2; j++) {
                file >> K[i][j];
            }
        }

        file.close();
    }

    MPI_Bcast(&t, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D1, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&D2, 1, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);

    if (rank != 0) {
        K.resize(D1, vector<int>(D2));
    }
    for (int i = 0; i < D1; i++) {
        MPI_Bcast(K[i].data(), D2, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    }

    int rowsPerProc = n / size;
    int row_padding = (D1 - 1) / 2; // notie ////////////////////////////////////
    int col_padding = (D2 - 1) / 2;
    int startRow = rank * rowsPerProc;
    int endRow = (rank == size - 1) ? n : startRow + rowsPerProc;

    vector<vector<int>> local_A(endRow - startRow + 2 * row_padding, vector<int>(m, 0));
    new_A.resize(endRow - startRow + 2 * row_padding, vector<int>(m, 0));

    if (rank == 0) {
        // Send
        for (int r = 1; r < size; r++) {
            int sendStartRow = r * rowsPerProc;
            int sendEndRow = (r == size - 1) ? n : (r + 1) * rowsPerProc;

            for (int i = sendStartRow; i < sendEndRow; i++) {
                MPI_Send(A[i].data(), m, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD);
            }
        }

        for (int i = 0; i < endRow - startRow; i++) {
            local_A[i + row_padding] = A[i];
        }
    } else {
        // Receive
        for (int i = 0; i < endRow - startRow; i++) {
            MPI_Recv(local_A[i + row_padding].data(), m, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }

    for (int time = 0; time < t; time++) {
        vector<int> uppersend, lowersend, upperrecv(row_padding * m), lowerrecv(row_padding * m);
        for (int i = row_padding; i < 2 * row_padding; i++) {
            for (int j = 0; j < m; j++) {
                uppersend.push_back(local_A[i][j]);
            }
        }

        for (int i = local_A.size() - row_padding-1; i < local_A.size(); i++) {
            for (int j = 0; j < m; j++) {
                lowersend.push_back(local_A[i][j]);
            }
        }

        MPI_Send(uppersend.data(), m * row_padding, MPI_LONG_LONG_INT, (rank - 1 + size) % size, 0, MPI_COMM_WORLD);
        MPI_Send(lowersend.data(), m * row_padding, MPI_LONG_LONG_INT, (rank + 1) % size, 1, MPI_COMM_WORLD);

        MPI_Recv(upperrecv.data(), m * row_padding, MPI_LONG_LONG_INT, (rank - 1 + size) % size, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        MPI_Recv(lowerrecv.data(), m * row_padding, MPI_LONG_LONG_INT, (rank + 1) % size, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        for (int i = 0; i < row_padding; i++) {
            for (int j = 0; j < m; j++) {
                local_A[i][j] = upperrecv[(i)*m + j];
            }
        }
        for (int i = 0; i < row_padding; i++) {
            for (int j = 0; j < m; j++) {
                local_A[local_A.size() - row_padding + i][j] = lowerrecv[i * m + j];
            }
        }

        for (int i = row_padding; i < local_A.size() - row_padding; i++) {
            for (int j = 0; j < m; j++) {
                int sum = 0;
                for (int di = -row_padding; di <= row_padding; di++) {
                    for (int dj = -col_padding; dj <= col_padding; dj++) {
                        int x = i + di;
                        int y = (j + dj + m) % m;
                        sum += K[di + row_padding][dj + col_padding] * local_A[x][y];
                    }
                }
                new_A[i][j] = sum / (D1 * D2);
            }
        }

        local_A = new_A;

        MPI_Barrier(MPI_COMM_WORLD);
    }

    // result
    if (rank == 0) {
        for (int r = 1; r < size; r++) {
            int recvStartRow = r * rowsPerProc;
            int recvEndRow = (r == size - 1) ? n : (r + 1) * rowsPerProc;
            for (int i = recvStartRow; i < recvEndRow; i++) {
                MPI_Recv(A[i].data(), m, MPI_LONG_LONG_INT, r, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            }
        }

        for (int i = 0; i < endRow - startRow; i++) {
            A[i] = local_A[i + row_padding];
        }
    } else {
        for (int i = 0; i < endRow - startRow; i++) {
            MPI_Send(local_A[i + row_padding].data(), m, MPI_LONG_LONG_INT, 0, 0, MPI_COMM_WORLD);
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < m; j++) {
                cout << A[i][j] << " ";
            }
        }
    }

    MPI_Finalize();
    return 0;
}
