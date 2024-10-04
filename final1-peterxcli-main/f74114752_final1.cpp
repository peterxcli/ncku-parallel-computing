#include<bits/stdc++.h>
using namespace std;

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

int main() {
    string input_file_name = "";
    cin >> input_file_name;
    
    ifstream file; 
    file.open(input_file_name);
    int n, k;
    file >> n >> k;
    vector<int> a(n);
    for (auto &x : a) file >> x;
    auto res = counting(a, k);
    for (auto &x : res) cout << x << " ";
}