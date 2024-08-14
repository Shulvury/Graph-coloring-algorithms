#include "Graph.h"
#include "inputGraph.h"
#include "reactcol.h"
#include "tabu.h"
#include "manipulateArrays.h"
#include "initializeColoring.h"
#include <iomanip>
#include <string.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <limits.h>
#include<random>

//This makes sure the compiler uses _strtoui64(x, y, z) with Microsoft Compilers, otherwise strtoull(x, y, z) is used
#ifdef _MSC_VER
#define strtoull(x, y, z) _strtoui64(x, y, z)
#endif

using namespace std;

unsigned long long numConfChecks;

std::random_device rd;
std::mt19937 mt(rd());

unsigned long long maxChecks = 100000000;

void logAndExitNow(string s) {
	//Writes message s to screen and log file and then exits the program
	ofstream resultsLog("resultsLog.log", ios::app);
	resultsLog << s;
	cout << s;
	resultsLog.close();
	exit(1);
}

void makeVectorAdjList(vector< vector<int> >& adjList, Graph& g)
{
	//Makes the adjacency list corresponding to G
	adjList.clear();
	for (int i = 0; i < g.n; i++) {
		adjList.push_back(vector<int>());
	}
	for (int i = 0; i < g.n; i++) {
		for (int j = 0; j < g.n; j++) {
			if (g[i][j] && i != j) {
				adjList[i].push_back(j);
			}
		}
	}
}

int SpanNewEvalue(Graph& g, vector<int>& col, vector<int>& T, vector< vector<int> >& adjList, int& k) {

	int val = 0;
	for (int i = 0; i < g.n; i++) {
		for (auto v : adjList[i]) {
			for (auto t : T) {
				numConfChecks++;
				if (col[v] - col[i] == t) {
					val++;
					break;
				}
			}
		}
	}
	
	return val;
}

vector< vector<int> > SpanCalculateMatrix(Graph& g, vector<int>& col, vector<int>& T, int& k, vector <vector<int> >& adjList) {
	
	vector< vector<int> > M;
	for (int i = 0; i < g.n; i++) {
		M.push_back(vector<int>(k, 0));
	}
	numConfChecks += g.n * k;

	//now set up conflict matrix M
	int Tsize = T.size();

	for (int i = 0; i < g.n; i++) {
		for (int c = 0; c < k; c++) {
			for (auto v : adjList[i]) {
				for (auto t : T) {
					numConfChecks++;
					if (c - col[v] == t) {
						M[i][c]++;
						break;
					}
				}
			}
		}
	}

	return M;
}


int SpanTabuCol(Graph& g, vector<int>& col, vector<int>& T, vector< vector<int> >& adjList, int& k, vector<vector<int>>& M, int& E) {

	uniform_int_distribution<int> r_t(0, 9);

	vector<vector<int>> tabuList;
	for (int i = 0; i < g.n; i++) {
		tabuList.push_back(vector<int>(k,-1));
	}

	numConfChecks += g.n * k;

	int bestGlobalVal = E;
	int currentIteration = 0;
	

	while (numConfChecks <= maxChecks) {

		vector<int> VertInConfl;
		for (int i = 0; i < g.n; i++) {
			if (M[i][col[i]] > 0) VertInConfl.push_back(i);
		}
		numConfChecks += g.n;

		int bestVal = INT_MAX;
		pair<int, int> bestMove = make_pair(-1, -1);
		int newval = 0;

		for (auto i : VertInConfl) {														//move all vertices in conflict
			for (int c = 0; c < k; c++) {													//into all colors
				if (c == col[i]) continue;													// except its own color
				newval = E - 2*M[i][col[i]] + 2*M[i][c];										// calculate new Energy state
				numConfChecks++;
				if (newval < bestVal) {														// if its better than current best
					numConfChecks += 2;
					if (tabuList[i][c] < currentIteration || newval < bestGlobalVal) {		// , not in tabu, or the very best globally
						bestMove.first = i;													// save it as best move
						bestMove.second = c;
						bestVal = newval;													// and value as current best value
						if (newval < bestGlobalVal) bestGlobalVal = newval;
					}
				}
			}
		}

		// now best move is found, time to execute it
		if (bestMove.first == -1) { //if all moves are tabu???, take random move
			uniform_int_distribution<int> rv(0, g.n - 1);
			uniform_int_distribution<int> rk(0, k - 1);
			bestMove.first = rv(mt);
			bestMove.second = rk(mt);
		}
		//first update M
		for (auto v : adjList[bestMove.first]) {
			numConfChecks += 2*T.size();
			for (auto t : T) {
				if (col[bestMove.first] - col[v] == t) {											//the vertex is moved from prev color
					M[bestMove.first][col[bestMove.first]]--;
					M[v][col[v]]--;
				}	
				if (bestMove.second - col[v] == t) {													//the vertex is moved to new color
					M[bestMove.first][bestMove.second]++;
					M[v][col[v]]++;
				}
			}
		}
		//update tabu tenure list
		int nc = 0;
		for (int i = 0; i < g.n; i++) {
			if (nc > M[i][col[i]]) nc = M[i][col[i]];
		}
		numConfChecks += g.n + 2;
		tabuList[bestMove.first][col[bestMove.first]] = currentIteration + 0.6 * nc + r_t(mt);
		//tabuList[bestMove.first][col[bestMove.first]] = currentIteration + r_t(mt);

		//now change color
		col[bestMove.first] = bestMove.second;

		if (bestVal == 0) return 0;

		currentIteration++;
	}

	return bestGlobalVal;
}

void SpanInitialDSatur(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, vector<int>& T) {

	int i, j, r;

	//Make a vector representing all the nodes
	vector<int> permutation(g.n);
	for (i = 0; i < g.n; i++)permutation[i] = i;
	//Randomly permute the nodes, and then arrange by increasing order of degree
	//(this allows more than 1 possible outcome from the sort procedure)
	for (i = permutation.size() - 1; i >= 0; i--) {
		r = rand() % (i + 1);
		swap(permutation[i], permutation[r]);
	}
	//Bubble sort is used here. This could be made more efficent
	for (i = (permutation.size() - 1); i >= 0; i--) {
		for (j = 1; j <= i; j++) {
			numConfChecks += 2;
			if (adjList[permutation[j - 1]].size() > adjList[permutation[j]].size()) {
				swap(permutation[j - 1], permutation[j]);
			}
		}
	}

	//We also have a vector to hold the saturation degrees of each node
	vector<int> satDeg(permutation.size(), 0);

	//Initialise candSol and colNode
	for (i = 0; i < colNode.size(); i++) colNode[i] = INT_MIN;

	//Colour the rightmost node first (it has the highest degree), and remove it from the permutation
	colNode[permutation.back()] = 0;
	int nodePos = permutation.back();
	permutation.pop_back();
	//..and update the saturation degree array
	satDeg.pop_back();
	for (i = 0; i < satDeg.size(); i++) {
		numConfChecks++;
		if (g[nodePos][permutation[i]]) {
			satDeg[i]++;
		}
	}

	//Now colour the remaining nodes.
	nodePos = 0;
	int maxSat;
	while (!permutation.empty()) {
		//choose the node to colour next (the rightmost node that has maximal satDegree)
		maxSat = INT_MIN;
		for (i = 0; i < satDeg.size(); i++) {
			if (satDeg[i] >= maxSat) {
				maxSat = satDeg[i];
				nodePos = i;
			}
		}
		//now choose which colour to assign to the node
		int clr = g.n * g.n;
		bool foundColor = true;
		for (int i = 0; i < clr; i++) {
			foundColor = true;
			for (auto v : adjList[nodePos]) {
				if (!foundColor) break;
				for (auto t : T) {
					if (colNode[v] - i == t && colNode[v] > -1) {
						foundColor = false;
					}
				}
			}
			if (foundColor) {
				colNode[nodePos] = i;
				break;
			}
		}


		//Finally, we remove the node from the permutation
		permutation.erase(permutation.begin() + nodePos);
		satDeg.erase(satDeg.begin() + nodePos);
	}

}

void SpanDSaturKcol(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, int k, vector<int>& T) {

	int i, j, r;
	bool foundColour;

	//Make a vector representing all the nodes
	vector<int> permutation(g.n);
	for (i = 0; i < g.n; i++)permutation[i] = i;
	//Randomly permute the nodes, and then arrange by increasing order of degree
	//(this allows more than 1 possible outcome from the sort procedure)
	for (i = permutation.size() - 1; i >= 0; i--) {
		r = rand() % (i + 1);
		swap(permutation[i], permutation[r]);
	}
	//Bubble sort is used here. This could be made more efficent
	for (i = (permutation.size() - 1); i >= 0; i--) {
		for (j = 1; j <= i; j++) {
			numConfChecks += 2;
			if (adjList[permutation[j - 1]].size() > adjList[permutation[j]].size()) {
				swap(permutation[j - 1], permutation[j]);
			}
		}
	}

	//We also have a vector to hold the saturation degrees of each node
	vector<int> satDeg(permutation.size(), 0);

	//Initialise candSol and colNode
	vector< vector<int> > candSol;
	for (int i = 0; i < k; i++) {
		candSol.push_back(vector<int>());
	}
	for (i = 0; i < colNode.size(); i++) colNode[i] = INT_MIN;

	//Colour the rightmost node first (it has the highest degree), and remove it from the permutation
	candSol[0].push_back(permutation.back());
	colNode[permutation.back()] = 0;
	permutation.pop_back();
	//..and update the saturation degree array
	satDeg.pop_back();
	for (i = 0; i < satDeg.size(); i++) {
		numConfChecks++;
		if (g[candSol[0][0]][permutation[i]]) {
			satDeg[i]++;
		}
	}

	//Now colour the remaining nodes.
	int nodePos = 0, maxSat;
	while (!permutation.empty()) {
		//choose the node to colour next (the rightmost node that has maximal satDegree)
		maxSat = INT_MIN;
		for (i = 0; i < satDeg.size(); i++) {
			if (satDeg[i] >= maxSat) {
				maxSat = satDeg[i];
				nodePos = i;
			}
		}
		//now choose which colour to assign to the node
		foundColour = false;
		int freq = -1;

		for (int c = 0; c < candSol.size(); c++) { //for all colors
			foundColour = true;
			for (int v = 0; v < adjList[nodePos].size(); v++) { //for all adjacency list
				if (!foundColour) break;
				for (int t = 0; t < T.size(); t++) {
					if (!foundColour) break;
					if (c - colNode[v] == T[t] && colNode[v] > -1) foundColour = false; //check if assigning C to V makes any illegal difference
				}
			}
			if (foundColour) {
				freq = c;
				break;
			}
		}
		if (foundColour) {
			candSol[freq].push_back(nodePos);
			colNode[nodePos] = freq;
		}

		//Finally, we remove the node from the permutation
		permutation.erase(permutation.begin() + nodePos);
		satDeg.erase(satDeg.begin() + nodePos);
	}

	uniform_int_distribution<int> uintd(0, k - 1);

	for (int i = 0; i < g.n; i++) {
		//assign random color to other vertices
		if (colNode[i] < 0) {
			colNode[i] = uintd(mt);
		}
	}
}








int main(int argc, char** argv)
{

	Graph g;
	int k, frequency = 0, increment = 0, verbose = 0, randomSeed = 1, tenure = 1, algorithm = 2, cost, duration, constructiveAlg = 1, targetCols = 1;


	//Program parameters

	inputDimacsGraph(g, "dsjc1000.9.col");

	//This variable keeps count of the number of times information about the instance is looked up 
	numConfChecks = 0;

	//Seed
	srand(randomSeed);


	//Make the adjacency list structure 
	int** neighbors = new int* [g.n];
	makeAdjList(neighbors, g);


	// span structures
	vector<vector<int>> adjList;
	makeVectorAdjList(adjList, g);

	vector<int> col(g.n, 0);


	//The solution is held in the following array
	int* coloring = new int[g.n];
	int* bestColouring = new int[g.n];

	int mode = 2;
	/* SELECT MODE
		1 -- SPAN TABU
		2 -- PLAIN GRAPH COL TABU
	*/

	//vector<int> T = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 14, 15};
	vector<int> T = { 0};

	//Now start the timer
	clock_t clockStart = clock();

	//Generate the initial value for k using greedy or dsatur algorithm
	if (mode == 1) {
		SpanInitialDSatur(col, g, adjList, T);
		k = *max_element(col.begin(), col.end()) + 1;
	}
	else if (mode == 2) {
		k = generateInitialK(g, constructiveAlg, bestColouring);
		//..and write the results to the output file
	}

	cout << "\tInitial coloring gave " << k << " colors. Starting tabu search.\n\n";

	//MAIN ALGORITHM
	k--;
	
	maxChecks = 1e9;
	

	if (mode == 1) {
		while (numConfChecks < maxChecks) {

			//Do the algorithm for this value of k, either until a slution is found, or maxChecks is exceeded
			SpanDSaturKcol(col, g, adjList, k, T);

			vector<vector<int>> M;
			M = SpanCalculateMatrix(g, col, T, k, adjList);
			
			int E = 0;
			E = SpanNewEvalue(g, col, T, adjList, k);
			cout << "\t inital cost is " << E << " while edges are " << g.nbEdges << endl;

			cost = SpanTabuCol(g, col, T, adjList, k, M, E);

			//Algorithm has finished at this k
			if (cost == 0) {

				//Copy the current solution as the best solution
				for (int i = 0; i < g.n; i++) bestColouring[i] = col[i];
				//Check if the target has been met
				if (k <= targetCols) {
					cout << "Reached desired number of colors" << k << endl;
					break;
				}
				cout << "\tReached " << k << " proper coloring using " << numConfChecks / 1000 << "k checks. Continuing\n";
				k--;
				continue;
			}
			cout << "  MaxChecks exceeded. Best found proper coloring has " << k + 1 << " colors. Terminating\n";
			cout << "\t\t Best found cost is " << cost << endl;
		}
	}
	else if (mode == 2) {
		while (numConfChecks < maxChecks) {

			//Initialise the solution array
			for (int i = 0; i < g.n; i++) coloring[i] = 0;

			//Do the algorithm for this value of k, either until a slution is found, or maxChecks is exceeded
			cost = tabu(g, coloring, k, maxChecks, tenure, verbose, frequency, increment, neighbors);

			//Algorithm has finished at this k
			if (cost == 0) {

				for (int i = 0; i < g.n; i++) {
					for (int j = i + 1; j < g.n; j++) {
						if (g[i][j] > 0 && coloring[i] == coloring[j]) cout << "Achtung!!!!\n";
					}
				}

				//Copy the current solution as the best solution
				for (int i = 0; i < g.n; i++) bestColouring[i] = coloring[i] - 1;
				//Check if the target has been met
				if (k <= targetCols) {
					cout << "Reached desired number of colors" << k << endl;
					break;
				}
				cout << "\tReached " << k << " proper coloring using " << numConfChecks / 1000 << "k checks. Continuing\n";
				duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
				cout << "Time taken: " << duration << "\n\n";
				k--;
				continue;
			}
			cout << "  MaxChecks exceeded. Best found proper coloring has " << k + 1 << " colors. Terminating\n";
		}
	}


	duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
	cout << "\tTook total of " << duration << " miliseconds\n";

	//Delete the arrays and end	
	delete[] coloring;
	delete[] bestColouring;
	for (int i = 0; i < g.n; i++) delete[] neighbors[i];
	delete[] neighbors;


	return 0;
}



