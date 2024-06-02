#include "Solver.hpp"

Solver::Solver(string fileName, int burninglength, int sample_size, int verbose) {
    this->burninglength = burninglength;
    this->readGraph(fileName);
    this->calculateMiddleNodes((verbose & (1 << 0)));
    this->findComponents();
    this->betweennessCentrality(fileName, sample_size, (verbose & (1 << 1)));
}

void Solver::readGraph(string graphFileName) {
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

void Solver::findComponents() {
    nodeMark = vector<int>(nodes, 0);
    componentsSize = vector<int>(nodes, 0);

    for (int i = 0; i < nodes; i++) {
        if (nodeMark[i] != 0)
            continue;
        
        vector<int> component;
        for (int j = 0; j < nodes; j++) {
            if (distanceMatrix[i][j] < USHRT_MAX) {
                nodeMark[j] = i + 1;
                component.push_back(j);
            }
        }
        components.push_back(component);
        for (int j = 0; j < nodes; j++) {
            if (distanceMatrix[i][j] < USHRT_MAX) {
                componentsSize[j] = component.size();
            }
        }
    }
}

void Solver::betweennessCentrality(string graphFileName, int sampleSize, int verbose) {   
    vector<double> centralityValues(nodes, 0.0), componenetsMaxCentrality(nodes, 0.0);

    graphFileName = "./Data/" + graphFileName;
    centralityValues = run_RAND1(graphFileName, min(sampleSize, nodes));
    for (int i = 0; i < nodes; i++)
        componenetsMaxCentrality[nodeMark[i]] = max(componenetsMaxCentrality[nodeMark[i]], centralityValues[i]);
    
    for (int i = 0; i < nodes; i++) {
        if (componenetsMaxCentrality[nodeMark[i]] == 0)
            centrality.push_back({1, i});
        else
            centrality.push_back({centralityValues[i] / componenetsMaxCentrality[nodeMark[i]], i});
    }
    
    sort(centrality.begin(), centrality.end(), greater<pair<double, int>>());

    if (verbose)
        for (auto nodeValues : centrality) {
            cout << nodeValues.first << " " << nodeValues.second << endl;
            if (nodeValues.first < 0.01)
                break;
        }
}

void Solver::BFSDistance(int start) {
    queue<int> BFSQueue;
    nodeMark = vector<int>(nodes, 0);
    vector<unsigned short int> height(nodes, USHRT_MAX);    
    
    BFSQueue.push(start);
    height[start] = 0;
    nodeMark[start] = 1;

    while (BFSQueue.size()) {
        int u = BFSQueue.front();
        BFSQueue.pop();
        for (int i = 0; i < adjList[u].size(); i++) {
            int v = adjList[u][i];
            if (nodeMark[v] == 0) {
                nodeMark[v] = 1;
                height[v] = height[u] + 1;
                BFSQueue.push(v);
            }
        }
    }
    distanceMatrix[start] = height;
}

void Solver::BFSMiddle(int start) {
    queue<int> BFSQueue;
    vector<int> seq;
    nodeMark = vector<int>(nodes, 0);
    vector<unsigned short int> middle(nodes, USHRT_MAX);    
    vector<unsigned short int> parent(nodes, 0);    
    
    BFSQueue.push(start);
    seq.push_back(start);
    nodeMark[start] = 1;

    while (BFSQueue.size()) {
        int u = BFSQueue.front();
        BFSQueue.pop();
        for (int i = 0; i < adjList[u].size(); i++) {
            int v = adjList[u][i];
            if (nodeMark[v] == 0) {
                parent[v] = u;
                nodeMark[v] = 1;
                BFSQueue.push(v);
                seq.push_back(v);
            }
        }
    }
    int i = 0, ptr = 0;
    while (i < seq.size()) {
        int u = seq[i]; int v = seq[ptr];
        if (distanceMatrix[start][u] % 2 == 1) {
            middle[u] = middle[parent[u]];
            i++;
        }
        else if (((int)distanceMatrix[start][v] == (int)distanceMatrix[v][u]) &&
                ((int)distanceMatrix[start][v] + (int)distanceMatrix[v][u] == (int)distanceMatrix[start][u])) {
            middle[u] = v;
            i++;
        } else {
            ptr++;
        }
    }
    middleMatrix[start] = middle;
}

void Solver::calculateMiddleNodes(int printToFile) {
    middleMatrix = vector<vector<unsigned short int>> (nodes, vector<unsigned short int>());
    distanceMatrix = vector<vector<unsigned short int>> (nodes, vector<unsigned short int>());

    for (int i = 0; i < nodes; i++)
        BFSDistance(i);
    for (int i = 0; i < nodes; i++)
        BFSMiddle(i);

    if (printToFile) {
        ofstream middleMatrixFile("middleMatrix.txt");
        for (auto vec : middleMatrix)
            for (int i = 0; i < vec.size(); i++)
                middleMatrixFile << vec[i] << " \n"[i == ((int)vec.size() - 1)];
        middleMatrixFile.close();

        ofstream distanceMatrixFile("distanceMatrix.txt");
        for (auto vec : distanceMatrix)
            for (int i = 0; i < vec.size(); i++)
                distanceMatrixFile << vec[i] << " \n"[i == ((int)vec.size() - 1)];
        distanceMatrixFile.close();
    }
}

vector<int> Solver::generateChromosome(int chromosomeSize, int minimumDistance, int verbose) {
    vector<int> chromosome(chromosomeSize, 0);
    mt19937 generator(randomDevice());
    for (int found = 0; found < chromosomeSize; found++) {
        vector<pair<double, int>> tempSigmoidCenrtality(nodes, {0, 0});
        double weightSum = 0;
        for (int i = 0; i < nodes; i++) {
            tempSigmoidCenrtality[i] = sigmoidCenrtality[i];
            for (int j = 0; j < found; j++)
                if (distanceMatrix[centrality[i].second][centrality[chromosome[j]].second] < minimumDistance)
                    tempSigmoidCenrtality[i].first = 0;
            weightSum += tempSigmoidCenrtality[i].first;
        }
        if (weightSum == 0) {
            minimumDistance--;
            found--;
            continue;
        }
        uniform_real_distribution<> nodeSelector(0, weightSum);
        double sampleNode = nodeSelector(generator);
        for (int i = 0; i < nodes; i++) {
            sampleNode -= tempSigmoidCenrtality[i].first;
            if (sampleNode < 0) {
                chromosome[found] = i;
                break;
            }
        }
    }
    vector<int> resultChromosome;
    for (auto node : chromosome)
        resultChromosome.push_back(centrality[node].second);
    
    if (verbose)
        for (auto node : resultChromosome)
            cout << node << " ";

    return resultChromosome;
}

void Solver::findNotBurnt(vector<int> chromosome, int verbose) {
    notBurnt.clear();
    nodeMark = vector<int>(nodes, nodes);
    nodeFillMark = vector<int>(nodes, nodes);
    
    for (int i = 0; i < chromosome.size(); i++)
        for (int node = 0; node < nodes; node++)
            nodeMark[node] = min(nodeMark[node], max(0, ((int)distanceMatrix[chromosome[i]][node] - (burninglength - 1 - i))));

    for (int node = 0; node < nodes; node++)
        if (nodeMark[node] != 0)
            notBurnt.push_back(node);
    
    if (verbose) {
        cout << notBurnt.size() << endl;
        for (auto node : notBurnt)
            cout << node << " ";
        cout << endl;
    }
}

vector<int> Solver::number2sequence(long long number, int sequencelength, int sequenceBase) {
    vector<int> sequence(sequencelength, 0);
    for (int i = 0; number > 0; i++) {
        sequence[i] = number % sequenceBase;
        number /= sequenceBase;
    }
    return sequence;
}

long long Solver::costFunction(vector<int> chromosome, int skipValue, int verbose) {
    findNotBurnt(chromosome, (verbose & (1 << 0)));

    if ((int)notBurnt.size() > skipValue)
        return INF;

    int notBurntSize = (int)notBurnt.size();
    int remaininglength = burninglength - (int)chromosome.size();

    if (notBurntSize == 0) {
        bestBurningSequence.clear();
        for (auto node : chromosome)
            bestBurningSequence.push_back(node);
        for (int i = 0; i < remaininglength; i++)
            bestBurningSequence.push_back(0);
        return 0;
    }

    long long statesNumber = 1;
    for (int i = 0; i < remaininglength; i++)
        statesNumber *= notBurntSize;
    
    long long minimumCost = INF;
    for (long long number = 0; number < statesNumber; number++) {
        vector<int> sequence = number2sequence(number, remaininglength, notBurntSize);

        long long cost = 0;
        for (auto node : notBurnt)
            nodeFillMark[node] = nodes;

        for (int i = 0; i < remaininglength; i++)
            for (auto node : notBurnt)
                nodeFillMark[node] = min(nodeFillMark[node], max(0, (int)distanceMatrix[node][notBurnt[sequence[i]]] - (remaininglength - 1 - i)));
        
        for (auto node : notBurnt)
            if (nodeFillMark[node])
                cost += min(nodeMark[node], nodeFillMark[node]) * min(nodeMark[node], nodeFillMark[node]);

        if ((verbose & (1 << 1)) && cost == 0) {
            for (auto nodeIndex : sequence)
                cout << notBurnt[nodeIndex] << " ";
            cout << endl;
        }
        minimumCost = min(minimumCost, cost);
        if (minimumCost == 0) {
            bestBurningSequence.clear();
            for (auto node : chromosome)
                bestBurningSequence.push_back(node);
            for (auto nodeIndex : sequence)
                bestBurningSequence.push_back(notBurnt[nodeIndex]);
            break;
        }
    }
    return minimumCost;
}

vector<int> Solver::crossover(vector<int> parentChromosome1, vector<int> parentChromosome2) {
    mt19937 generator(randomDevice());
    uniform_int_distribution<> methodSelector(0, 2);
    vector<int> childChromosome;
    
    for (int i = 0; i < (int)parentChromosome1.size(); i++) {
        switch (methodSelector(generator)) {
            case 0:
                if (middleMatrix[parentChromosome1[i]][parentChromosome2[i]] < nodes) {
                    childChromosome.push_back(middleMatrix[parentChromosome1[i]][parentChromosome2[i]]);
                    break;
                }
            case 1:
                childChromosome.push_back(parentChromosome1[i]);
                break;
            case 2:
                childChromosome.push_back(parentChromosome2[i]);
                break;
        }
    }
    return childChromosome;
}

vector<int> Solver::mutate(vector<int> chromosome, double mutateProbabilty) {
    mt19937 generator(randomDevice());
    uniform_real_distribution<> mutateSelector(0, 1);
    uniform_int_distribution<>  methodSelector(0, 1);

    for (int i = 0; i < (int)chromosome.size(); i++) {
        if (mutateSelector(generator) <= mutateProbabilty) {
            switch (methodSelector(generator)) {
                case 0: {
                    uniform_real_distribution<> nodeSelector(0, sigmoidCenrtalitySum);
                    double sampleNode = nodeSelector(generator);
                    for (int j = 0; j < nodes; j++) {
                        sampleNode -= sigmoidCenrtality[j].first;
                        if (sampleNode < 0) {
                            chromosome[i] = centrality[j].second;
                            break;
                        }
                    }
                }
                break;
                case 1: {
                    int neighborSize = adjList[chromosome[i]].size();
                    uniform_int_distribution<> neighborSelector(0, neighborSize - 1);
                    chromosome[i] = adjList[chromosome[i]][neighborSelector(generator)];
                }
                break;
            }
        }
    }
    return chromosome;
}

int Solver::rouletteWheelSelection(int generation, int topPopulation, double fitnessSum) {
    double maxWeight = 1.0 / ((double) population[generation][0].first + 1);
    mt19937 generator(randomDevice());
    uniform_int_distribution<> chromosomeSelector(0, population[generation].size() - 1);
    uniform_real_distribution<> chromosomeAcceptor(0, maxWeight);
    while (true) {
        int index = chromosomeSelector(generator);
        double weight = 1.0 / ((double) population[generation][index].first + 1);
        if (chromosomeAcceptor(generator) <= weight) {
            return index;
        }
    }
}

vector<int> Solver::solve(int chromosomeSize, int minimumDistance, 
                          int skipValue, int maxGenerations, int topPopulation,
                          int crossoverPopulation, double mutateProbabilty,
                          double alpha, double beta) {
    if (chromosomeSize == -1)
        chromosomeSize = burninglength - 3;

    if (minimumDistance == -1)
        minimumDistance = burninglength - 3;

    sigmoidCenrtality = vector<pair<double, int>>(nodes, {0, 0});
    sigmoidCenrtalitySum = 0;
    for (int i = 0; i < nodes; i++) {
        sigmoidCenrtality[i].first = 1.0 / (exp(-1 * beta * (centrality[i].first - alpha)) + 1);
        sigmoidCenrtality[i].second = centrality[i].second;
        sigmoidCenrtalitySum += sigmoidCenrtality[i].first;
    }

    for (int i = 0; i < topPopulation; i++) {
        cout << '*'; fflush(stdout);
        vector<int> chromosome = generateChromosome(chromosomeSize, minimumDistance);
        int fitness = costFunction(chromosome, skipValue);
        population[0].push_back({fitness, chromosome});
        if (fitness == 0)
            break;
    }
    cout << endl;

    for (int generation = 0; generation < maxGenerations; generation++) {
        sort(population[generation].begin(), population[generation].end());
        population[generation].erase(unique(population[generation].begin(), population[generation].end()), population[generation].end());
        
        if (population[generation][0].first == 0)
            return bestBurningSequence;

        if (population[generation][0].first == INF)
            skipValue += 10;

        while ((int)population[generation].size() < topPopulation) {
            vector<int> chromosome = generateChromosome(chromosomeSize, minimumDistance);
            int fitness = costFunction(chromosome, skipValue);
            population[generation].push_back({fitness, chromosome});            
        }
        sort(population[generation].begin(), population[generation].end());
        
        double fitnessSum = 0;
        for (int i = 0; i < topPopulation; i++)
            fitnessSum += population[generation][i].first;

        cout << "#" << generation + 1 << ": " << endl;
        for (int i = 0; i < 10; i++)
            cout << population[generation][i].first << " ";
        cout << endl;

        for (int i = 0; i < topPopulation; i++)
            population[generation + 1].push_back(population[generation][i]);
        
        for (int i = 0; i < crossoverPopulation; i++) {
            int chromosomeParent1 = 0, chromosomeParent2 = 0;
            while (chromosomeParent1 == chromosomeParent2) {
                chromosomeParent1 = rouletteWheelSelection(generation, topPopulation, fitnessSum);
                chromosomeParent2 = rouletteWheelSelection(generation, topPopulation, fitnessSum);
            }
            vector<int> chromosome = crossover(population[generation][chromosomeParent1].second, 
                                               population[generation][chromosomeParent2].second);
            int fitness = costFunction(chromosome, skipValue);
            population[generation + 1].push_back({fitness, chromosome});
        }
        for (int i = 0; i < topPopulation + crossoverPopulation; i++) {
            vector<int> chromosome = mutate(population[generation + 1][i].second, mutateProbabilty);  
            int fitness = costFunction(chromosome, skipValue);
            population[generation + 1].push_back({fitness, chromosome});
        }
    }
    return bestBurningSequence;
}
