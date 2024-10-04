#pragma GCC optimize("O2")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("fast-math")
#pragma GCC optimize("no-stack-protector")
#include <bits/stdc++.h>
#include <pthread.h>

using namespace std;

// #define int long long
vector<vector<int> > A_t, A_t1, K;
int n, m, D1, D2, row_padding, col_padding;

int thread_num = 0;
int CHUNK_SIZE = 69;

void *transform(void *arg) {
    int start_row = *(int *)arg;

    for (int row = start_row; row <= start_row + CHUNK_SIZE && row < n; ++row) {
        for (int j = 0; j < m; ++j) {
         int sum = 0;
            for (int di = -row_padding; di <= row_padding; di++) {
                for (int dj = -col_padding; dj <= col_padding; dj++) {
                    int wrapped_i = (row + di + n) % n;
                    int wrapped_j = (j + dj + m) % m;
                    sum += K[di + row_padding][dj + col_padding] * A_t[wrapped_i][wrapped_j];
                }
            }
            A_t1[row][j] = sum / (D1 * D2);
        }
    }

    return NULL;
}

int32_t main(int argc, char *argv[]) {
    thread_num = atoi(argv[1]);
    cin.tie(0);
    ios::sync_with_stdio(false);
    cout.tie(0);
    string filename;
    cin >> filename;
    ifstream fin(filename);
    int t;
    fin >> t;
    fin >> n >> m;
    CHUNK_SIZE = n / thread_num + 1;
    A_t.resize(n, vector<int>(m, 0));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m; j++) {
            fin >> A_t[i][j];
        }
    }
    fin >> D1 >> D2;
    row_padding = D1 / 2, col_padding = D2 / 2;
    K.resize(D1, vector<int>(D2, 0));
    for (int i = 0; i < D1; i++) {
        for (int j = 0; j < D2; j++) {
            fin >> K[i][j];
        }
    }

    A_t1.resize(n, vector<int>(m, 0));
    pthread_t threads[thread_num];
    int start_indices[thread_num];
    for (int time = 0; time < t; time++) {
        int current_thread = 0;
        for (int i = 0; i < n; i += CHUNK_SIZE) {
            start_indices[current_thread] = i;
            pthread_create(&threads[current_thread], NULL, transform, (void *)&start_indices[current_thread]);

            current_thread++;
            if (current_thread == thread_num || i + CHUNK_SIZE >= n) {
                for (int j = 0; j < current_thread; ++j) {
                    pthread_join(threads[j], NULL);
                }
                current_thread = 0;
            }
        }
        for (int i = 0; i < n; i++) for (int j = 0; j < m; j++) A_t[i][j] = A_t1[i][j];
    }
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j)
            cout << A_t1[i][j] << " ";
    }
}