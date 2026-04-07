#include <iostream>
#include <vector>
#include <algorithm>
#include <fstream>

using namespace std;

// Struct to represent a monomial like x_a * x_b * x_c
struct Mono {
    int a, b, c;
    bool operator==(const Mono& o) const { return a==o.a && b==o.b && c==o.c; }
    bool operator<(const Mono& o) const {
        if (a != o.a) return a < o.a;
        if (b != o.b) return b < o.b;
        return c < o.c;
    }
};

// Recursive function to search combinations and check Lex-Min
void find_orbits(int start, int depth, vector<int>& comb, const vector<vector<int>>& action_table, ofstream& out) {
    if (depth == comb.size()) {
        bool is_canonical = true;
        for (int pi = 0; pi < 720; ++pi) {
            vector<int> mapped(comb.size());
            for (size_t i = 0; i < comb.size(); ++i) {
                mapped[i] = action_table[pi][comb[i]];
            }
            sort(mapped.begin(), mapped.end());
            
            // Lexicographical comparison
            bool smaller = false;
            for (size_t i = 0; i < comb.size(); ++i) {
                if (mapped[i] < comb[i]) { smaller = true; break; }
                if (mapped[i] > comb[i]) { break; }
            }
            if (smaller) {
                is_canonical = false;
                break;
            }
        }
        if (is_canonical) {
            for (size_t i = 0; i < comb.size(); ++i) {
                out << comb[i] << (i + 1 == comb.size() ? "" : ",");
            }
            out << "\n";
        }
        return;
    }
    for (int i = start; i < 56; ++i) {
        comb[depth] = i;
        find_orbits(i + 1, depth + 1, comb, action_table, out);
    }
}

int main() {
    int K_VAL = 6; // Change this for k=6, k=7, etc.
    string filename = "orbit_reps_k" + to_string(K_VAL) + ".txt";
    
    // 1. Generate 56 Monomials (degree 3 in 6 variables)
    vector<Mono> monos;
    for(int i=0; i<6; ++i)
        for(int j=i; j<6; ++j)
            for(int k=j; k<6; ++k)
                monos.push_back({i,j,k});
                
    // 2. Generate 720 permutations of S6
    vector<vector<int>> perms;
    vector<int> p = {0,1,2,3,4,5};
    do { perms.push_back(p); } while(next_permutation(p.begin(), p.end()));
    
    // 3. Build Action Table
    vector<vector<int>> action_table(720, vector<int>(56));
    for(int pi=0; pi<720; ++pi) {
        for(int mi=0; mi<56; ++mi) {
            int arr[3] = {perms[pi][monos[mi].a], perms[pi][monos[mi].b], perms[pi][monos[mi].c]};
            sort(arr, arr+3);
            Mono mapped = {arr[0], arr[1], arr[2]};
            auto it = lower_bound(monos.begin(), monos.end(), mapped);
            action_table[pi][mi] = distance(monos.begin(), it);
        }
    }
    
    // 4. Run Lex-Min Search
    cout << "Finding orbits for k=" << K_VAL << "..." << endl;
    ofstream out(filename);
    vector<int> comb(K_VAL);
    find_orbits(0, 0, comb, action_table, out);
    out.close();
    
    cout << "Done! Saved to " << filename << endl;
    return 0;
}