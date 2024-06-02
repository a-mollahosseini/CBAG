#include "Checker.hpp"

Checker::Checker(string fileName) {
    this->readGraph(fileName);
}

void Checker::readGraph(string graphFileName) {
    graphFileName = "./Data/" + graphFileName;
    cout << graphFileName << endl;
    ifstream graphFile(graphFileName);

    string graphName; graphFile >> graphName; 
    graphFile >> nodes >> edges;
    int oneIndexed; graphFile >> oneIndexed;
    int weighted; graphFile >> weighted;

    vector<pair<int, int>> edgeList;
    for (int i = 0; i < edges; i++) {
        int u, v; double w;
        graphFile >> u >> v;
        if (oneIndexed) {
            u--; v--;
        }
        if (weighted == 1)
            graphFile >> w;

        if (u > v)
            swap(u, v);
        edgeList.push_back({u, v});
    }
    graphFile.close();

    sort(edgeList.begin(), edgeList.end());
    edgeList.erase(unique(edgeList.begin(), edgeList.end()), edgeList.end());

    edges = (int)edgeList.size();
    adjList = vector<vector<int>>(nodes, vector<int>());

    for (int i = 0; i < edges; i++) {
        int u = edgeList[i].first;
        int v = edgeList[i].second;
        adjList[u].push_back(v);
        adjList[v].push_back(u);
    }
    cout << graphName << " " << nodes << " " << edges << endl;
}

bool Checker::check(vector<int> solution, int verbose) {
    int len = (int)solution.size();
    vector<bool> nodeMark(nodes, 0);
    queue<pair<int, int>> BFSQueue;
    int count = 1, burnIndex = 0;

    BFSQueue.push({0, solution[0]});
    nodeMark[solution[0]] = 1;
    
    while (BFSQueue.size()) {
        pair<int, int> u = BFSQueue.front(); 
        BFSQueue.pop();
        
        for (auto v : adjList[u.second]) {
            int depth = u.first + 1;
            if (nodeMark[v] == 0 && depth < len) {
                BFSQueue.push({depth, v});
                nodeMark[v] = 1;
                count++;
                if (depth > burnIndex) {
                    burnIndex = depth;
                    if (nodeMark[solution[burnIndex]] == 0) {
                        BFSQueue.push({depth, solution[burnIndex]});
                        nodeMark[solution[burnIndex]] = 1;
                        count++;
                    }
                }
            }
        }
        if ((int)BFSQueue.size() == 0 && (burnIndex + 1) < len) {
            burnIndex++;
            int depth = burnIndex;
            BFSQueue.push({depth, solution[burnIndex]});
            if (nodeMark[solution[burnIndex]] == 0) {
                nodeMark[solution[burnIndex]] = 1;
                count++;
            }
        }
    }
    if (count == nodes)
        return true;
    else {
        if (verbose) {
            for (int node = 0; node < nodes; node++)
                if (!nodeMark[node])
                    cout << node << " ";
            cout << endl;
        }
        return false;
    }
}
