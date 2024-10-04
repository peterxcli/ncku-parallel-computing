#include <mpi.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define INF 1000000000
#define TAG_DATA 0

short w[2500000000];

void readInput(char *filename, long long *n) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        perror("Failed to open file");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }
    fscanf(f, "%lld", n);

    long long x, y;
    short c;
    while (fscanf(f, "%lld %lld %hd", &x, &y, &c) != EOF) {
        w[y * (*n) + x] = c + 1;
    }
    fclose(f);
}

void distributeData(long long n, int rank, int world_size, short **loc_w, int *loc_begin, int *loc_end) {
    long long displs = n;
    for (int i = 0; i < world_size; ++i) {
        int sendcount = ((n - 1) / world_size + (i < (n - 1) % world_size)) * n;
        if (rank == 0 && i != 0) {
            MPI_Send(w + displs, sendcount, MPI_SHORT, i, TAG_DATA, MPI_COMM_WORLD);
        }

        if (rank == i) {
            *loc_w = (short *)malloc(sendcount * sizeof(short));
            if (!*loc_w) {
                perror("Failed to allocate memory");
                MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
            }
            *loc_begin = displs / n;
            *loc_end = *loc_begin + sendcount / n;

            if (i != 0) {
                MPI_Recv(*loc_w, sendcount, MPI_SHORT, 0, TAG_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            } else {
                memcpy(*loc_w, w + displs, sendcount * sizeof(short));
            }
        }
        displs += sendcount;
    }
}

void initializeDistances(long long n, int loc_n, short *loc_w, int *dist, bool *vis) {
    for (int i = 0; i < loc_n; ++i) {
        dist[i] = loc_w[i * n] == 0 ? INF : loc_w[i * n] - 1;
        vis[i] = false;
    }
}

void dijkstra(long long n, int rank, int loc_begin, int loc_end, short *loc_w, int *dist, bool *vis) {
    int min_pair[2], global_min_pair[2];
    while (true) {
        // Find vertex with minimum distance
        int ui = -1;
        for (int i = 0; i < loc_end - loc_begin; ++i) {
            if (!vis[i] && (ui == -1 || dist[i] < dist[ui])) {
                ui = i;
            }
        }

        // Break if all vertices are visited
        if (ui == -1) break;

        // Reduce global minimum
        min_pair[0] = dist[ui];
        min_pair[1] = ui + loc_begin;
        MPI_Allreduce(min_pair, global_min_pair, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        int u = global_min_pair[1], min_dist = global_min_pair[0];

        // Mark u as visited if it is in this process
        if (u >= loc_begin && u < loc_end) {
            vis[u - loc_begin] = true;
        }

        // Update distance
        for (int i = 0; i < loc_end - loc_begin; ++i) {
            if (!vis[i] && loc_w[i * n + u] != 0 &&
                dist[i] > min_dist + loc_w[i * n + u]) {
                dist[i] = min_dist + loc_w[i * n + u] - 1;
            }
        }
    }

    // Keep reducing until all processes have finished
    while (global_min_pair[1] != -1) {
        min_pair[0] = INF;
        min_pair[1] = -1;
        MPI_Allreduce(min_pair, global_min_pair, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
    }
}

void outputResults(int rank, int world_size, int loc_n, int *dist) {
    if (rank == 0) printf("0 ");
    for (int i = 0; i < world_size; ++i) {
        if (rank == i) {
            for (int j = 0; j < loc_n; ++j) {
                printf("%d ", dist[j]);
            }
        }
        MPI_Barrier(MPI_COMM_WORLD);
    }
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);

    int rank, world_size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    long long n;
    if (rank == 0) {
        char filename[100];
        scanf("%s", filename);
        readInput(filename, &n);
    }

    MPI_Bcast(&n, 1, MPI_LONG_LONG, 0, MPI_COMM_WORLD);

    short *loc_w;
    int loc_begin, loc_end;
    distributeData(n, rank, world_size, &loc_w, &loc_begin, &loc_end);

    int loc_n = loc_end - loc_begin;
    int *dist = (int *)malloc(loc_n * sizeof(int));
    bool *vis = (bool *)malloc(loc_n * sizeof(bool));
    initializeDistances(n, loc_n, loc_w, dist, vis);

    dijkstra(n, rank, loc_begin, loc_end, loc_w, dist, vis);

    outputResults(rank, world_size, loc_n, dist);
}
