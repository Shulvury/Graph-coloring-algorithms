#ifndef INITIALIZECOLORING_INCLUDED
#define INITIALIZECOLORING_INCLUDED

#include "Graph.h"
#include <vector>

int generateInitialK(Graph& g, int alg, int* bestColouring);
void initializeColoring(Graph& g, int* c, int k);
void initializeColoringForTabu(Graph& g, int* c, int k);

void assignAColourDSatur(bool& foundColour, std::vector< std::vector<int> >& candSol, std::vector<int>& permutation, int nodePos, std::vector<int>& satDeg, Graph& g, std::vector<int>& colNode, std::vector< std::vector<int> >& adjList);


#endif
