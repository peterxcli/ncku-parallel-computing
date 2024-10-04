#pragma GCC optimize("Ofast")
#pragma GCC optimize("inline")
#pragma GCC optimize("-fgcse")
#pragma GCC optimize("-fgcse-lm")
#pragma GCC optimize("-fipa-sra")
#pragma GCC optimize("-ftree-pre")
#pragma GCC optimize("-ftree-vrp")
#pragma GCC optimize("-fpeephole2")
#pragma GCC optimize("-ffast-math")
#pragma GCC optimize("-fsched-spec")
#pragma GCC optimize("unroll-loops")
#pragma GCC optimize("-falign-jumps")
#pragma GCC optimize("-falign-loops")
#pragma GCC optimize("-falign-labels")
#pragma GCC optimize("-fdevirtualize")
#pragma GCC optimize("-fcaller-saves")
#pragma GCC optimize("-fcrossjumping")
#pragma GCC optimize("-fthread-jumps")
#pragma GCC optimize("-funroll-loops")
#pragma GCC optimize("-fwhole-program")
#pragma GCC optimize("-freorder-blocks")
#pragma GCC optimize("-fschedule-insns")
#pragma GCC optimize("inline-functions")
#pragma GCC optimize("-ftree-tail-merge")
#pragma GCC optimize("-fschedule-insns2")
#pragma GCC optimize("-fstrict-aliasing")
#pragma GCC optimize("-fstrict-overflow")
#pragma GCC optimize("-falign-functions")
#pragma GCC optimize("-fcse-skip-blocks")
#pragma GCC optimize("-fcse-follow-jumps")
#pragma GCC optimize("-fsched-interblock")
#pragma GCC optimize("-fpartial-inlining")
#pragma GCC optimize("no-stack-protector")
#pragma GCC optimize("-freorder-functions")
#pragma GCC optimize("-findirect-inlining")
#pragma GCC optimize("-fhoist-adjacent-loads")
#pragma GCC optimize("-frerun-cse-after-loop")
#pragma GCC optimize("inline-small-functions")
#pragma GCC optimize("-finline-small-functions")
#pragma GCC optimize("-ftree-switch-conversion")
#pragma GCC optimize("-foptimize-sibling-calls")
#pragma GCC optimize("-fexpensive-optimizations")
#pragma GCC optimize("-funsafe-loop-optimizations")
#pragma GCC optimize("inline-functions-called-once")
#pragma GCC optimize("-fdelete-null-pointer-checks")
#pragma GCC optimize("O2,unroll-loops")
#pragma GCC target("avx2,bmi,bmi2,lzcnt,popcnt")
#include <bits/stdc++.h>
#include <mpi.h>
using namespace std;

struct Point {
    int x, y;
    int num;
    bool operator<(const Point &rhs_p) const {
        if (x != rhs_p.x) {
            return x < rhs_p.x;
        } else {
            return y < rhs_p.y;
        }
    }
    pair<int, int> operator-(const Point &rhs_p) const {
        return {x - rhs_p.x, y - rhs_p.y};
    }
};

template <typename T>
T cross(const pair<T, T> &a, const pair<T, T> &b) {
    return a.first * b.second - a.second * b.first;
}

template <typename T>
vector<T> getConvexHull(vector<T> &pnts) {
    int n = pnts.size();
    sort(pnts.begin(), pnts.end());
    vector<T> hull;
    for (int i = 0; i < 2; i++) {
        // i = 0: find lower hull
        // i = 1: find upper hull
        int t = hull.size();
        for (T pnt : pnts) {
            while (hull.size() - t >= 2 && cross(hull.back() - hull[hull.size() - 2], pnt - hull[hull.size() - 2]) >= 0)
                hull.pop_back();
            hull.emplace_back(pnt);
        }
        hull.pop_back();
        reverse(pnts.begin(), pnts.end());
    }
    return hull;
}
int mod(int a, int b) {
    return (a % b + b) % b;
}

vector<Point> merge_convex_hull(const vector<Point> &left_hull, const vector<Point> &right_hull) {
    int left_hull_rightmost_index = max_element(left_hull.begin(), left_hull.end()) - left_hull.begin();
    int right_hull_leftmost_index = min_element(right_hull.begin(), right_hull.end()) - right_hull.begin();

    // Find upper tangent
    int left_hull_upper_index = left_hull_rightmost_index, right_hull_upper_index = right_hull_leftmost_index;
    bool left_done = false, right_done = false;
    while (!left_done || !right_done) {
        // fix right, move left to upper
        auto old_new_tangent_cross =
            cross(
                right_hull[right_hull_upper_index] - left_hull[left_hull_upper_index],
                left_hull[mod(left_hull_upper_index - 1, left_hull.size())] - left_hull[left_hull_upper_index]);
        if (old_new_tangent_cross > 0) {
            left_hull_upper_index = mod(left_hull_upper_index - 1, left_hull.size());
            left_done = false;
        } else { // the left point of the upper tangent is found
            left_done = true;
        }

        // fix left, move right to upper
        old_new_tangent_cross =
            cross(
                left_hull[left_hull_upper_index] - right_hull[right_hull_upper_index],
                right_hull[mod(right_hull_upper_index + 1, right_hull.size())] - right_hull[right_hull_upper_index]);
        if (old_new_tangent_cross < 0) {
            right_hull_upper_index = mod(right_hull_upper_index + 1, right_hull.size());
            right_done = false;
        } else {
            right_done = true;
        }
    }

    // Find lower tangent
    int left_hull_lower_index = left_hull_rightmost_index, right_hull_lower_index = right_hull_leftmost_index;
    left_done = false, right_done = false;
    while (!left_done || !right_done) {
        // fix right, move left to lower
        auto old_new_tangent_cross =
            cross(
                right_hull[right_hull_lower_index] - left_hull[left_hull_lower_index],
                left_hull[mod(left_hull_lower_index + 1, left_hull.size())] - left_hull[left_hull_lower_index]);

        if (old_new_tangent_cross <= 0) {
            left_hull_lower_index = mod(left_hull_lower_index + 1, left_hull.size());
            left_done = false;
        } else { // the left point of the upper tangent is found
            left_done = true;
        }

        // fix left, move right to lower
        old_new_tangent_cross =
            cross(
                left_hull[left_hull_lower_index] - right_hull[right_hull_lower_index],
                right_hull[mod(right_hull_lower_index - 1, right_hull.size())] - right_hull[right_hull_lower_index]);
        if (old_new_tangent_cross >= 0) {
            right_hull_lower_index = mod(right_hull_lower_index - 1, right_hull.size());
            right_done = false;
        } else {
            right_done = true;
        }
    }

    // merge convex hull
    int left_hull_leftmost_index = min_element(left_hull.begin(), left_hull.end()) - left_hull.begin();
    vector<Point> convex_hull;
    for (int i = left_hull_leftmost_index; i != left_hull_upper_index; i = mod(i + 1, left_hull.size())) {
        convex_hull.emplace_back(left_hull[i]);
    }
    convex_hull.emplace_back(left_hull[left_hull_upper_index]);

    for (int i = right_hull_upper_index; i != right_hull_lower_index; i = mod(i + 1, right_hull.size())) {
        convex_hull.emplace_back(right_hull[i]);
    }
    convex_hull.emplace_back(right_hull[right_hull_lower_index]);

    for (int i = left_hull_lower_index; i != left_hull_leftmost_index; i = mod(i + 1, left_hull.size())) {
        convex_hull.emplace_back(left_hull[i]);
    }
    return convex_hull;
}

int32_t main(int32_t argc, char **argv) {
    int32_t rank, size;
    MPI_Init(&argc, &argv);               // Initialize the MPI execution environment
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the calling process in the communicator
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the size of the group associated with a communicatoroi

    MPI_Datatype mpi_point_type;
    MPI_Type_contiguous(3, MPI_INT, &mpi_point_type);
    MPI_Type_commit(&mpi_point_type);

    int n;
    string filename;
    vector<Point> points;
    vector<Point> local_points;
    
    if (rank == 0) {
        cin >> filename;
        ifstream fin(filename);
        fin >> n;
        points.resize(n);
        for (int i = 0; i < n; i++) {
            points[i].num = i+1;
            fin >> points[i].x >> points[i].y;
        }
        sort(points.begin(), points.end());
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

    local_points.resize(sendcounts[rank]);
    MPI_Scatterv(points.data(), sendcounts.data(), displs.data(),
               mpi_point_type, local_points.data(), sendcounts[rank],
               mpi_point_type, 0, MPI_COMM_WORLD);

    vector<Point> local_convex_hull = getConvexHull(local_points);

    for (int step = 1; step < size; step *= 2) {
        if (rank % (2 * step) == 0) {
            // Even ranked process: receive convex hull from right neighbor and merge
            if (rank + step < size) {
                int received_size;
                MPI_Recv(&received_size, 1, MPI_INT, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                vector<Point> received_points(received_size);
                MPI_Recv(&received_points[0], received_size * sizeof(Point), MPI_BYTE, rank + step, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                
                local_convex_hull = merge_convex_hull(local_convex_hull, received_points);
            }
        } else {
            // Odd ranked process: send convex hull to left neighbor
            int sending_size = local_convex_hull.size();
            MPI_Send(&sending_size, 1, MPI_INT, rank - step, 0, MPI_COMM_WORLD);
            MPI_Send(&local_convex_hull[0], sending_size * sizeof(Point), MPI_BYTE, rank - step, 0, MPI_COMM_WORLD);
            break;
        }
    }

    if (rank == 0) {
        // Final convex hull is in rank 0
        for (const auto& point : local_convex_hull) {
            cout << point.num << " ";
        }
        // cout << endl;
    }

    MPI_Finalize(); // Terminates MPI execution environment
    return 0;
}
