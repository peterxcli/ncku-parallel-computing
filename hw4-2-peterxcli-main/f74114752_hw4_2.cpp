#include <bits/stdc++.h>
#include <pthread.h>

using namespace std;
#define int long long

int thread_num;
vector<int> p;
int **m;
int n;
int L;

typedef struct {
    int start;
    int end;
} ThreadParam;

void *computeMatrixChain(void *param) {
    ThreadParam *threadParam = (ThreadParam *)param;
    int start = threadParam->start, end = threadParam->end;

    for (int i = start; i < min(end, n - L); i++) {
        int j = i + L;
        m[i][j] = LLONG_MAX;
        for (int k = i; k < j; k++) {
            int q = m[i][k] + m[k + 1][j] + p[i] * p[k + 1] * p[j + 1];
            if (q < m[i][j])
                m[i][j] = q;
        }
    }

    pthread_exit(NULL);
}

int MatrixChainOrder(int n) {
    m = new int *[n + 1];
    for (int i = 0; i <= n; i++) {
        m[i] = new int[n + 1];
        m[i][i] = 0;
    }

    for (L = 1; L < n; L++) {
        pthread_t threads[thread_num];
        ThreadParam params[thread_num];
        int chunk = (n - L) / thread_num;
        for (int i = 0; i < thread_num; i++) {
            params[i].start = i * chunk;
            params[i].end = (i == thread_num - 1) ? n - L : (i + 1) * chunk;
            pthread_create(&threads[i], NULL, computeMatrixChain, (void *)&params[i]);
        }

        for (int i = 0; i < thread_num; i++) {
            pthread_join(threads[i], NULL);
        }
    }

    return m[0][n - 1];
}

int32_t main(int32_t argc, char *argv[]) {
    thread_num = atoi(argv[1]);
    string filename;
    cin >> filename;
    ifstream fin(filename);
    fin >> n;
    p.resize(n + 1);
    for (int i = 1; i <= n; i++)
        fin >> p[i];
    p[0] = 1;

    cout << MatrixChainOrder(n);

    for (int i = 0; i <= n; i++) {
        delete[] m[i];
    }
    delete[] m;

    return 0;
}
