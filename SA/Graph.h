#ifndef GraphIncluded
#define GraphIncluded

class Graph {
public:

	Graph();
	Graph(int n);
	~Graph();
	//  Graph(char * file);

	void resize(int n);
	void assignWeight(int w);

	int* matrix;
	int n;        // number of nodes
	int nbEdges;  // number of edges

	int* operator[](int index);
};

#endif
