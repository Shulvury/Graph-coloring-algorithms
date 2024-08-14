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
#include <vector>
#include <limits.h>
#include <random>
#include <set>
#include "SA.h"

//This makes sure the compiler uses _strtoui64(x, y, z) with Microsoft Compilers, otherwise strtoull(x, y, z) is used
#ifdef _MSC_VER
#define strtoull(x, y, z) _strtoui64(x, y, z)
#endif

using namespace std;

unsigned long long numConfChecks;
unsigned long long maxChecks = 100000000;
clock_t clockStart;

//uniform_real_distribution<double> unif(0.0, 1.0);
//random_device rd;
//mt19937 mt(rd());

double alpha = 0.99; // cooling rate for temperature
double t_min;




int main(int argc, char** argv){

	Graph g;
	int k, frequency = 0, increment = 0, verbose = 0, randomSeed = 2, tenure = 1, algorithm = 2, cost, duration, constructiveAlg = 1, targetCols = 2;

	inputDimacsGraph(g, "C2000.5.col.txt");


	/* MODES OF OPERATION
		1 == SPAN (incomplete)
		2 == normal SA for graph col
		3 == FAP for matrix contraint
		4 == Modified Penalty function
		5 == Original pemalty function
	*/
	int mode = 2;


	numConfChecks = 0;

	maxChecks = 1e9;
	
	//Seed
	srand(randomSeed);


	int* first_col = new int[g.n];
	for (int i = 0; i < g.n; i++) {
		first_col[i] = 0;
	}

	vector<int> col(g.n);
	vector< vector<int> > adjList;
	makeAdjList(adjList, g);
	numConfChecks += g.nbEdges;

	vector<int> T = {0,1,2,3,4};

	

	//Now start the timer
	clockStart = clock();


	if (mode == 1) {
		delete[] first_col;
		//SpanInitialDSatur(col, g, adjList, T);
		SpanInitialRandomCol(g, col, T, adjList);
		k = *max_element(col.begin(), col.end()) + 1;

		bool colorTargeFailed = false;

		cout << "DSatur coloring gave " << k << " colors.\n";

		while (numConfChecks < maxChecks && k + 1 > targetCols) {

			int E = 0; // "energy" of system that we want to minimize
			double temp = 100.0; // starting temperature;
			int loop = 100;		 // tries per one temp before cooling
			alpha = 0.95;		 // cooling rate
			t_min = 1.0;		 // minimal temperature, stop when reached

			maxChecks = 500000000;

			k--;

			
			SpanDSaturKcol(col, g, adjList, k, T); // new improper coloring with one less color
			SpanCalculateEnergy(g, col, adjList, E, T);

			temp = E * 5;

			cout << "starting energy is " << E << endl;

			SpanSimulatedAnnealing(g, col, k, adjList, temp, E, loop, T);

			cout << "end energy is " << E << endl;
			SpanCalculateEnergy(g, col, adjList, E, T);
			cout << "checking span calculate energy function, E = " << E << endl << endl;

			//for (int i = 0; i < g.n; i++) {
			//	cout << i << " vert --- col " << col[i] << endl;
			//}
			//int confl_count = 0;
			//for (int i = 0; i < g.n; i++) {
			//	for (auto v : adjList[i]) {
			//		if (col[i] - col[v] == 0) {
			//			cout << "vertices " << v << " and " << i << " have a same color!\n";
			//			confl_count++;
			//		}
			//	}
			//}

			SpanCalculateEnergy(g, col, adjList, E, T);
			cout << "checking span calculate energy function, E = " << E << endl << endl;
			//cout << "\t conflict number is " << confl_count/2 << endl;

			if (E == 0) {
				cout << "\tSuccessfully acquired proper " << k << "-coloring. Total checks: " << numConfChecks / 1000 << "k\n";
				colorTargeFailed = false;
			}
			else if (numConfChecks < maxChecks) {
				cout << "\tBest proper coloring is still " << k + 1 << endl;
				cout << "\tTrying again.\n\n";
				colorTargeFailed = true;
				k++;
			}
			else {
				cout << "\tBest proper coloring found has " << k + 1 << endl;
				cout << "\Terminating.\n";
			}
		}
	}
	else if (mode == 2) {

		//Generate the initial value for k using dsatur algorithm
		k = generateInitialK(g, 1, first_col);
		for (int i = 0; i < g.n; i++) {
			col[i] = first_col[i];
		}

		// show result
		cout << "DSatur coloring gave " << k << " colors.\n";
		delete[] first_col;

		//k = 105;
		//MAIN ALGORITHM

		while (numConfChecks < maxChecks && k + 1 > targetCols) {

			int E = 0;				 // "energy" of system that we want to minimize
			double temp = 6.0;	 // starting temperature;
			int loop = 4*g.n;			 // tries per one temp before cooling
			alpha = 0.95;			 // cooling rate
			t_min = 0.000001;			 // minimal temperature, stop when reached


			k--;

			//conflict matrix 
			vector<vector<int>> M;
			for (int i = 0; i < g.n; i++) {
				M.push_back(vector<int>(k));
			}

			//DSaturKcol(col, g, adjList, k); // new improper coloring with one less color
			randomizeLastColor(g, col, k);
			ConflictMatrix(g, col, k, M, adjList); // new matrix for k coloring
			vector<vector<int>> ColClass;
			for (int i = 0; i < k; i++) {
				ColClass.push_back({});
			}
			for (int i = 0; i < g.n; i++) {
				ColClass[col[i]].push_back(i);
			}
			
			//CalculateEnergy(g, col, adjList, E);
			for (int i = 0; i < g.n; i++) {
				E += M[i][col[i]];
			}
			E = E / 2;

			cout << "Staring with energy " << E << endl;
			SimulatedAnnealing(g, col, k, adjList, M, ColClass, temp, E, loop);
			cout << "Ending current run with energy " << E << endl;

			//CalculateEnergy(g, col, adjList, E);
			//cout << "\t Further checking shows, E = " << E << "\n";

			if (E == 0) {
				cout << "\tSuccessfully acquired proper " << k << "-coloring. Total checks: " << numConfChecks / 1000 << "k\n";
				for (int i = 0; i < g.n; i++) {
					for (int j = 0; j < g.n; j++) {
						if (g[i][j] && col[i] == col[j]) cout << "Houston we have a problem!\n";
					}
				}
				duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
				cout << "Time taken: " << duration << "\n\n";
				if (k == 2) return 0;
			}
			else if (numConfChecks < maxChecks) {
				cout << "\tBest proper coloring is still " << k + 1 << endl;
				cout << "\tTrying again.\n\n";
				k++;
			}
			else {
				cout << "\tBest proper coloring found has " << k + 1 << endl;
				cout << "\Terminating.\n";
			}
		}
	}
	else if (mode == 3) {
		delete[] first_col;
		
		//g.assignWeight(5);

		FAP_Intial_DSatur(col, g, adjList);
		k = *max_element(col.begin(), col.end()) + 1;

		cout << "DSatur coloring gave " << k << " colors.\n";

		while (numConfChecks < maxChecks && k + 1 > targetCols) {

			int E = 0; // "energy" of system that we want to minimize
			double temp = 100.0; // starting temperature;
			int loop = 100;		 // tries per one temp before cooling
			alpha = 0.95;		 // cooling rate
			t_min = 0.1;		 // minimal temperature, stop when reached

			maxChecks = 1e9;

			k--;


			FAP_K_col_Dsatur(col, g, adjList, k); // new improper coloring with one less color

			/*
			for (int i = 0; i < k + 10; i++) {
				cout << "Color " << i << ": ";
				for (int v = 0; v < g.n; v++) {
					if (col[v] == i) cout << v << ' ';
				}
				cout << "\n";
			}
			cout << "thats all\n\n";
			*/

			/*//conflict matrix 
			vector<vector<int>> M;
			for (int i = 0; i < g.n; i++) {
				M.push_back(vector<int>(k));
			}

			FAP_calcMatrix(g, col, adjList, M, k);
			E = FAP_energy(g, col, adjList, M, k);*/

			temp = FAP_Seek_T0(g, col, k, adjList);
			cout << "T_0 was set to " << temp << endl;
			E = FAP_InitialEnergy(g, col, adjList);
			cout << "starting energy is " << E << endl;
			if (E == 0) {
				cout << "\tSuccessfully acquired proper " << k << "-coloring. Total checks: " << numConfChecks / 1000 << "k\n";
				continue;
			}

			FAP_SimulatedAnnealing(g, col, k, adjList, temp, E, loop);

			cout << "end energy is " << E << endl;
			//cout << "lets check energy via calc function: " << FAP_energy(g, col, adjList, M, k) << endl;
			cout << "lets check energy via calc function: " << FAP_InitialEnergy(g, col, adjList) << endl;

			/*for (int v = 0; v < g.n; v++) {
				for (auto u : adjList[v]) {
					if (g[v][u] > abs(col[v] - col[u])) {
						cout << v << "--" << u << endl;
					}
				}
			}*/
			/*
			for (int i = 0; i < g.n; i++) {
				for (int c = 0; c < k; c++) {
					cout << M[i][c] << ' ';
				}
				cout << endl;
			}
			cout << "\n\n";
			*/
			if (E == 0) {
				cout << "\tSuccessfully acquired proper " << k << "-coloring. Total checks: " << numConfChecks / 1000 << "k\n";
			}
			else if (numConfChecks < maxChecks) {
				cout << "\tBest proper coloring is still " << k + 1 << endl;
				cout << "\tTrying again.\n\n";
				k++;
			}
			else {
				cout << "\tBest proper coloring found has " << k + 1 << endl;
				cout << "\Terminating.\n";
			}
		}
	}
	else if (mode == 4) {

		//Generate the initial value for k using dsatur algorithm
		k = generateInitialK(g, 1, first_col);
		for (int i = 0; i < g.n; i++) {
			col[i] = first_col[i];
		}

		// show result
		cout << "DSatur coloring gave " << k << " colors.";
		delete[] first_col;

		k--;
		//k = 105;
		//MAIN ALGORITHM

		while (numConfChecks < maxChecks) {

			int E = 0;				 // "energy" of system that we want to minimize
			double temp = 100.0;	 // starting temperature;
			int loop = 32*g.n;			 // tries per one temp before cooling
			alpha = 0.99;			 // cooling rate
			t_min = -1;			 // minimal temperature, stop when reached

			maxChecks = 1e9;

			//k--;

			//DSaturKcol(col, g, adjList, k); // new improper coloring with one less color
			randomizeLastColor(g, col, k);
			cout << '\n';

			vector<vector<int>> ColClass;
			for (int i = 0; i < k; i++) {
				ColClass.push_back({});
			}
			for (int i = 0; i < g.n; i++) {
				ColClass[col[i]].push_back(i);
			}
			

			E = PenaltySimulatedAnnealing(g, col, adjList, ColClass, temp, loop, k);
			if (E == 0) {
				for (int i = 0; i < g.n; i++) {
					for (auto v : adjList[i]) {
						if (col[i] == col[v]) cout << "FAILURE!!!\n";
					}
				}
				cout << "\t\t\t\tSucceded in aquarining " << k << " coloring!\n";
				k--;
			}

			duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
			cout << "Run ended. Time taken: " << duration << " || temp: " << temp << "\n";
			//t_min = temp;
			cout << "Currently at " << numConfChecks / 1000 << "k checks\n\n";

		}
	}
	else if (mode == 5) {

		//Generate the initial value for k using dsatur algorithm
		k = generateInitialK(g, 2, first_col);
		for (int i = 0; i < g.n; i++) {
			col[i] = first_col[i];
		}

		//k = g.nbEdges / g.n;
		//RandomCol(g, col, k);
	
		// show result
		cout << "DSatur coloring gave " << k << " colors.";
		delete[] first_col;

		
		//k = 105;
		//MAIN ALGORITHM

		int E = 0;				 // "energy" of system that we want to minimize
		double temp = 10.0;	 // starting temperature;
		int loop = 32 * g.n;			 // tries per one temp before cooling
		alpha = 0.95;			 // cooling rate
		t_min = -1;			 // minimal temperature, stop when reached

		maxChecks = 1e10;

		//DSaturKcol(col, g, adjList, k-1); // new improper coloring with one less color
		randomizeLastColor(g, col, k-1);
		cout << '\n';

		vector<vector<int>> ColClass;
		for (int i = 0; i < k; i++) {
			ColClass.push_back({});
		}
		for (int i = 0; i < g.n; i++) {
			ColClass[col[i]].push_back(i);
		}


		E = UnmodifiedPenaltySimulatedAnnealing(g, col, adjList, ColClass, temp, loop, k);

		duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
		cout << "Run ended. Time taken: " << duration << " || temp: " << temp << "\n";
	}


	duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);


	cout << "\n\tTook total of " << duration << " miliseconds.\n";

	return 0;
}



