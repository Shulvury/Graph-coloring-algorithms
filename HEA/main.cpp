#include "Graph.h"
#include "inputGraph.h"
#include "tabu.h"
#include "makesolution.h"
#include "xover.h"
#include "diversity.h"
#include <string.h>
#include <iostream>
#include <stdlib.h>
#include <fstream>
#include <stdio.h>
#include <time.h>
#include <limits.h>
#include <vector>
#include <iomanip>
#include <random>

//This makes sure the compiler uses _strtoui64(x, y, z) with Microsoft Compilers, otherwise strtoull(x, y, z) is used
#ifdef _MSC_VER
#define strtoull(x, y, z) _strtoui64(x, y, z)
#endif

using namespace std;

random_device rd;
mt19937 mt(rd());



void logAndExitNow(string s) {
	//Writes message s to screen and log file and then exits the program
	ofstream resultsLog("resultsLog.log", ios::app);
	resultsLog << s;
	cout << s;
	resultsLog.close();
	exit(1);
}

int ASS = INT_MIN;
bool solIsOptimal(vector<int>& sol, Graph& g, int k);
void makeAdjList(int** neighbors, Graph& g);
void replace(vector<vector<int> >& population, vector<int>& parents, vector<int>& osp, vector<int>& popCosts, Graph& g, int oCost);

unsigned long long numConfChecks;


int main(int argc, char** argv) {

	Graph g;
	uniform_int_distribution<int> rS(0, 1e7);
	int i, k, popSize = 10, maxIterations = 16, verbose = 0, randomSeed = 2, constructiveAlg = 1, targetCols = 2, xOverType = 1;
	bool solFound = false, doKempeMutation = false, measuringDiversity = false;
	unsigned long long maxChecks = 1e9;
	vector<int> parents;

	//This variable keeps count of the number of times information about the instance is looked up 
	numConfChecks = 0;

	randomSeed = rS(mt);

	inputDimacsGraph(g, "dsjc125.1.col");


	//Set the number of parents in each crossover and decide if the Kempe mutation is going to be used
	if (xOverType == 3) parents.resize(4);
	else parents.resize(2);
	if (xOverType == 2) doKempeMutation = true;

	//set tabucol limit
	maxIterations = maxIterations * g.n;
	if (targetCols < 2 || targetCols > g.n) targetCols = 2;

	//Now set up some output files

	//Make the adjacency list structure 
	int** neighbors = new int* [g.n];
	makeAdjList(neighbors, g);


	//Seed and start timer
	clock_t clockStart = clock();
	int duration;
	srand(randomSeed);
	numConfChecks = 0;

	//Data structures used for population and offspring
	vector<vector<int> > population(popSize, vector<int>(g.n));
	vector<int> popCosts(popSize);
	vector<int> osp(g.n), bestColouring(g.n);

	//Generate the initial value for k using greedy or dsatur algorithm
	k = generateInitialK(g, constructiveAlg, bestColouring);
	//..and write the results to the output file
	

	cout << "  Initial coloring gave " << k << " colors. Starting HEA\n";

	//MAIN ALGORITHM
	k--;
	while (numConfChecks < maxChecks && k + 1 > targetCols) {
		solFound = false;

		//First build the population
		for (i = 0; i < popSize; i++) {
			//Build a solution using modified DSatur algorithm
			makeInitSolution(g, population[i], k, verbose);
			//Check to see whether this solution is alrerady optimal or if the cutoff point has been reached. If so, we end
			if (solIsOptimal(population[i], g, k)) {
				solFound = true;
				for (int j = 0; j < g.n; j++)osp[j] = population[i][j];
				break;
			}
			if (numConfChecks >= maxChecks) {
				for (int j = 0; j < g.n; j++)osp[j] = population[i][j];
				break;
			}
			//Improve each solution via tabu search and record their costs
			popCosts[i] = tabu(g, population[i], k, maxIterations, 0, neighbors);
			//Check to see whether this solution is now optimal or if the cuttoff point is reached. If so, we end
			if (popCosts[i] == 0) {
				solFound = true;
				for (int j = 0; j < g.n; j++)osp[j] = population[i][j];
				break;
			}
			if (numConfChecks >= maxChecks) {
				for (int j = 0; j < g.n; j++)osp[j] = population[i][j];
				break;
			}
		}

		//Now evolve the population
		int rIts = 0, oCost = 1, best = INT_MAX;
		while (numConfChecks < maxChecks && !solFound) {

			//Choose parents and perform crossover to produce a new offspring 		
			doCrossover(xOverType, osp, parents, g, k, population);

			//Improve the offspring via tabu search and record its cost
			oCost = tabu(g, osp, k, maxIterations, 0, neighbors);

			//Write osp over weaker parent and update popCosts
			replace(population, parents, osp, popCosts, g, oCost);

			rIts++;

			if (oCost < best) best = oCost;
			if (oCost == 0) solFound = true;
		}

		//Algorithm has finished at this k	
		duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
		if (solFound) {
			//Copy the current solution as the best solution
			for (int i = 0; i < g.n; i++) bestColouring[i] = osp[i] - 1;

			cout << "\tSuccessfully found proper " << k << "-coloring. Used " << numConfChecks / 1000 << "k operations.\n\t Starting new iteration\n";
			duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
			cout << "Time spent: " << duration << " miliseconds.\n";
		}
		else {
			cout << "\t Failed. Best found proper coloring is " << k + 1 << " colors. Terminating\n";
		}

		k--;
	}


	duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
	cout << "\n\tTook total of " << duration << " miliseconds.\n";

	for (int i = 0; i < g.n; i++) delete[] neighbors[i];
	delete[] neighbors;

	return(0);
}

//*********************************************************************
inline
bool solIsOptimal(vector<int>& sol, Graph& g, int k)
{
	int i, j;
	for (i = 0; i < (g.n) - 1; i++) {
		for (j = i + 1; j < g.n; j++) {
			if (sol[i] == sol[j] && g[i][j])
				return(false);
		}
	}
	//If we are here then we have established a solution with k or fewer colours
	return(true);
}
//*********************************************************************
void makeAdjList(int** neighbors, Graph& g)
{
	//Makes the adjacency list corresponding to G
	for (int i = 0; i < g.n; i++) {
		neighbors[i] = new int[g.n + 1];
		neighbors[i][0] = 0;
	}
	for (int i = 0; i < g.n; i++) {
		for (int j = 0; j < g.n; j++) {
			if (g[i][j] && i != j) {
				neighbors[i][++neighbors[i][0]] = j;
			}
		}
	}
}

//*********************************************************************
void replace(vector<vector<int> >& population, vector<int>& parents, vector<int>& osp, vector<int>& popCosts, Graph& g, int oCost)
{
	//Go through the parents and ID the worst one
	int toDie = -1, i, max = INT_MIN;
	for (i = 0; i < parents.size(); i++) {
		if (popCosts[parents[i]] > max) {
			max = popCosts[parents[i]];
			toDie = parents[i];
		}
	}
	//Copy osp over the parent selected toDie
	population[toDie] = osp;
	popCosts[toDie] = oCost;
}
