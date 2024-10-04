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
    vector<vector<int>> A, K;
    string filename;

    if (rank == 0) {
        cin >> filename;
        ifstream file(filename);
        file >> t;
        file >> n >> m;
        A.resize(n, vector<int>(m));
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) file >> A[i][j];
        file >> D1 >> D2;
        K.resize(D1, vector<int>(D2));
        for (int i = 0; i < D1; i++) for (int j = 0; j < D2; j++) file >> K[i][j];
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

    if (rank != 0) {
        A.resize(n, vector<int>(m));
    }
    for (int i = 0; i < n; i++) {
        MPI_Bcast(A[i].data(), m, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
    }

    int chunk = n / (size-1);
    int row_padding = (D1 - 1) / 2;
    int col_padding = (D2 - 1) / 2;
    vector<int> sendcounts(size, chunk);
    vector<int> displs(size, 0);
    for (int i = 1; i < size; i++) {
        if (i == size - 1) {
            sendcounts[i] +=  n % (size-1);
        }
        displs[i] = (i - 1) * chunk;
    }
    int startRow = displs[rank];
    int endRow = startRow + sendcounts[rank];

    vector<vector<int>> local_A(n, vector<int>(m, 0));

    for (int time = 0; time < t; time++) {
        if (rank == 0) {
            for (int i = 1; i < size; i++) {
                int dis = displs[i];
                for (int j = dis; j < dis + sendcounts[i]; j++) {
                    MPI_Recv(A[j].data(), m, MPI_LONG_LONG_INT, i, j-displs[i] + 2, MPI_COMM_WORLD, (MPI_Status *)1);
                }
            }
        }
        else if (rank != 0) {
            for (int i = startRow; i < endRow; i++) {
                for (int j = 0; j < m; j++) {
                    int sum = 0;
                    for (int di = -row_padding; di <= row_padding; di++) {
                        for (int dj = -col_padding; dj <= col_padding; dj++) {
                            int x = (i + di + n) % n;
                            int y = (j + dj + m) % m;
                            sum += K[di + row_padding][dj + col_padding] * A[x][y];
                        }
                    }
                    local_A[i][j] = sum / (D1 * D2);
                }
            }

            for (int i = startRow; i < endRow; i++) {
                MPI_Send(local_A[i].data(), m, MPI_LONG_LONG_INT, 0, (i - startRow + 2), MPI_COMM_WORLD);
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
        for (auto &a_i : A) {
            MPI_Bcast(a_i.data(), m, MPI_LONG_LONG_INT, 0, MPI_COMM_WORLD);
        }
    }

    if (rank == 0) {
        for (const auto &row: A) {
            for (const auto ele: row) {
                cout << ele << " ";
            }
        }
    }

    MPI_Finalize();
    return 0;
}
