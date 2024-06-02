#include <bits/stdc++.h>

using namespace std;

class Checker {
    private:
        int nodes, edges;
        vector<vector<int>> adjList; 
    private:
        void readGraph(string graphFileName); 
    public:
        Checker(string fileName);
        bool check(vector<int> soloution, int verbose = 0);
};
