#include <algorithm>
#include <cfloat>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

using namespace std;

typedef struct {
    int x, y;
} Point;

bool operator<(const Point &a, const Point &b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
}

bool operator==(const Point &a, const Point &b) {
    return a.x == b.x && a.y == b.y;
}

int orientation(const Point &a, const Point &b, const Point &c) {
    int val = (b.y - a.y) * (c.x - b.x) - (b.x - a.x) * (c.y - b.y);

    if (val == 0)
        return 0;
    return val > 0 ? 1 : -1;
}

double dist(const Point &a, const Point &b) {
    double d = sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
    return floor(d * 10000) / 10000;
}

double mst(const vector<Point> points) {
    int n = points.size();
    vector<double> key(n, DBL_MAX);
    vector<bool> mst_set(n, false);
    vector<int> parent(n, -1);

    key[0] = 0;
    for (int i = 0; i < n - 1; ++i) {
        int u = -1;
        for (int j = 0; j < n; ++j) {
            if (!mst_set[j] && (u == -1 || key[j] < key[u])) {
                u = j;
            }
        }

        mst_set[u] = true;
        for (int v = 0; v < n; ++v) {
            if (!mst_set[v] && dist(points[u], points[v]) < key[v]) {
                parent[v] = u;
                key[v] = dist(points[u], points[v]);
            }
        }
    }

    double ans = 0;
    for (int i = 1; i < n; ++i) {
        ans += dist(points[i], points[parent[i]]);
    }
    return ans;
}

int main(int argc, char *argv[]) {
    vector<Point> convex_hull;
    vector<Point> inner_points;
    int n;
    string filename;
    cin >> filename;

    ifstream fin(filename);
    fin >> n;

    vector<Point> points(n);
    for (int i = 0; i < n; ++i) {
        fin >> points[i].x >> points[i].y;
    }

    sort(points.begin(), points.end());

    // Find convex hull
    for (int i = 0; i < 2; ++i) {
        int sz = convex_hull.size();
        for (const auto &point : points) {
            while (convex_hull.size() - sz >= 2 &&
                   orientation(convex_hull[convex_hull.size() - 2],
                               convex_hull[convex_hull.size() - 1], point) <= 0) {
                convex_hull.pop_back();
            }
            convex_hull.push_back(point);
        }

        convex_hull.pop_back();
        reverse(points.begin(), points.end());
    }

    // Find inner points
    for (const auto &point : points) {
        bool is_inner = true;
        for (const auto &hull_point : convex_hull) {
            if (point == hull_point) {
                is_inner = false;
                break;
            }
        }
        if (is_inner)
            inner_points.push_back(point);
    }

    // Find minimum spanning tree
    double ans = DBL_MAX;
    #pragma omp parallel for num_threads(8) reduction(min:ans)
    for (int i = 0; i < (1 << inner_points.size()); i += 1) {
        vector<Point> points(convex_hull);
        for (int j = 0; j < inner_points.size(); ++j) {
            if (i & (1 << j)) {
                points.push_back(inner_points[j]);
            }
        }
        ans = min(ans, mst(points));
    }

    // Reduce
    // double global_ans;
    // MPI_Reduce(&ans, &global_ans, 1, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
    printf("%.4lf", ans);

    // MPI_Finalize();
    return 0;
}