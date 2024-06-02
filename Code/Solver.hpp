#define BOOST_ALLOW_DEPRECATED_HEADERS
#define INF                             1000
#define MAX_NODES                       100000
#define MAX_GA_ITERATIONS               10000

#include <bits/stdc++.h>
#include "BetweennessApprox/RAND1.hpp"

using namespace std;

class Solver {
    private:
        int nodes, edges;
        int burninglength;
        double sigmoidCenrtalitySum;
        vector<int> nodeMark;
        vector<int> nodeFillMark;
        vector<vector<int>> adjList;
        vector<pair<double, int>> centrality;
        vector<pair<double, int>> sigmoidCenrtality;
        vector<vector<unsigned short int>> middleMatrix;
        vector<vector<unsigned short int>> distanceMatrix;
        vector<vector<int>> components;
        vector<int> componentsSize;
        vector<int> notBurnt;
        vector<pair<int, vector<int>>> population[MAX_GA_ITERATIONS];
        random_device randomDevice;
        vector<int> bestBurningSequence;
    private:
        void readGraph(string graphFileName);
        void findComponents();
        void betweennessCentrality(string graphFileName, int sampleSize, int verbose = 0);
        void BFSDistance(int start);
        void BFSMiddle(int start);
        void calculateMiddleNodes(int printToFile = 0);
        vector<int> generateChromosome(int chromosomeSize, int minimumDistance = 1, int verbose = 0);
        void findNotBurnt(vector<int> chromosome, int verbose = 0);
        vector<int> number2sequence(long long number, int sequencelength, int sequenceBase);
        long long costFunction(vector<int> chromosome, int skipValue, int verbose = 0);
        int rouletteWheelSelection(int generation, int topPopulation, double fitnessSum);
        vector<int> crossover(vector<int> parentChromosome1, vector<int> parentChromosome2);
        vector<int> mutate(vector<int> chromosome, double mutateProbabilty);
    public:
        Solver(string fileName, int burninglength, int sampleSize = 1000, int verbose = 0);
        vector<int> solve(int chromosomeSize = -1, int minimumDistance = -1, 
                          int skipValue = 20, int maxGenerations = 500, int topPopulation = 300,
                          int crossoverPopulation = 500, double mutateProbabilty = 0.1, 
                          double alpha = 0.05, double beta = 200);
};
