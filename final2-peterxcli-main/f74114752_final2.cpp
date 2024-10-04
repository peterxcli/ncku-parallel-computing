#include<bits/stdc++.h>
using namespace std;

int n;
vector<vector<int> > G;

vector<int> counting(const vector<int> &arr, int k) {
    vector<int> count(k+1), ret(arr.size());
    #pragma omp parallel for
    for (int i = 0; i < k; i++) {
        count[i] = 0;
    }
    #pragma omp parallel for
    for (int i = 0; i < arr.size(); i++) {
        count[arr[i]]++;
    }
    for (int i = 1; i < count.size(); i++) {
        count[i] = count[i] + count[i-1];
    }
    #pragma omp parallel for
    for (int i = arr.size()-1; i >= 0; i--) {
        ret[count[arr[i]]-1] = arr[i];
        count[arr[i]]--;
    }
    return ret;
}

vector<int> dijkstra() {
    vector<int> d(n, 10000000);
    priority_queue<pair<int32_t, int32_t>, vector<pair<int32_t, int32_t>>, greater<pair<int32_t, int32_t>>> pq;
    d[0] = 0;
    pq.push({0, 0});
    while (!pq.empty()) {
        auto [dist, u] = pq.top();
        pq.pop();
        #pragma omp parallel for num_threads(8)
        for (int v = 0; v < n; v++) {
            if (u == v) continue;
            if (d[u] + G[u][v] < d[v]) {
                d[v] = d[u] + G[u][v];
                #pragma omp critical
                {
                    pq.push({d[v], v});
                }
            }
        }
    }
    return d;
}

int main() {
    string input_file_name = "";
    cin >> input_file_name;
    
    ifstream file; 
    file.open(input_file_name);
    file >> n;
    G.resize(n, vector<int>(n));
    // vector<vector<int> > G(n, vector<int>(n));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file >> G[i][j];
        }
    }
    auto res = dijkstra(); 
    for (auto &x : res) cout << x << " ";
}