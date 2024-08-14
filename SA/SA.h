#pragma once

#include "Graph.h"
#include <vector>

using namespace std;

void makeAdjList(std::vector< std::vector<int> >& adjList, Graph& g);

void DSaturKcol(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, int k);

int SimanEnergy(Graph& g, vector<int>& col, vector<vector<int>>& ColClass, int k, vector<int>& E);

int SimanDeltaEnergy(Graph& g, vector<int>& col, vector<vector<int>>& ColClass, int node, int newcol, vector<int>& E);

void SpanInitialDSatur(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, vector<int>& T);

void SpanDSaturKcol(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, int k, vector<int>& T);

void SimulatedAnnealing(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList, vector< vector<int> >& M,
						vector<vector<int>>& ColClass, double& temp, int& E, int& loop);

void SpanSimulatedAnnealing(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList, double& temp, int& E, int& loop, vector<int>& T);

void ConflictMatrix(Graph& g, vector<int>& col, int k, vector< vector<int> >& M, vector< vector<int> >& adjList);

void CalculateEnergy(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int& E);

void SpanCalculateEnergy(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int& E, vector<int>& T);

void SpanInitialRandomCol(Graph& g, vector<int>& col, vector<int>& T, vector< vector<int>>& adjList);

void FAP_calcMatrix(Graph& g, vector<int>& col, vector< vector<int> >& adjList, vector<vector<int>>& M, int k);

int FAP_energy(Graph& g, vector<int>& col, vector< vector<int> >& adjList, vector<vector<int>>& M, int k);

void FAP_Execute_Move(Graph& g, vector<int>& col, vector< vector<int> >& adjList, vector<vector<int>>& M,
	int node, int newcolor, int& E, int k);

int FAP_EnergyDifference(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int node, int newcol);

void FAP_SimulatedAnnealing(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList,
	double& temp, int& E, int& loop);

void FAP_K_col_Dsatur(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, int k);

void FAP_Intial_DSatur(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList);

int FAP_InitialEnergy(Graph& g, vector<int>& col, vector<vector<int>>& adjList);

double FAP_Seek_T0(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList);

int PenaltySimulatedAnnealing(Graph& g, vector<int>& col, vector<vector<int>>& adjList, vector<vector<int>> ColClass,
	double& temp, int loop, int& k);

void randomizeLastColor(Graph& g, vector<int>& col, int k);

int UnmodifiedPenaltySimulatedAnnealing(Graph& g, vector<int>& col, vector<vector<int>>& adjList, vector<vector<int>> ColClass,
	double& temp, int loop, int& k);

void RandomCol(Graph& g, vector<int>& col, int k);
