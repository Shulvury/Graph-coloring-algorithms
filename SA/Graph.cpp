#include "Graph.h"
#include <fstream>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <random>

using namespace std;

Graph::Graph() {
	matrix = NULL;
	n = 0;
	nbEdges = 0;
}

Graph::Graph(int m) {
	matrix = NULL;
	resize(m);
}

int* Graph::operator[](int index) {
	if (index < 0 || index >= this->n) {
		cerr << "First node index out of range: " << index << "\n";
		matrix[-1] = 0; //Make it crash.
	}
	return matrix + this->n * index;
}

void Graph::resize(int m) {
	if (matrix != NULL) {
		delete[] matrix;
	}
	if (m > 0) {
		n = m;
		nbEdges = 0;
		matrix = new int[m * m];
		for (int i = 0; i < m * m; i++) {
			matrix[i] = 0;
		}
	}
}

void Graph::assignWeight(int w) {

	random_device randev;
	mt19937 merstw(randev());
	uniform_int_distribution<int> weight(1, w);

	if (matrix == NULL) return;
	for (int i = 0; i < n; i++) {
		for (int j = i + 1; j < n; j++) {
			if (i == j) continue;
			if (matrix[i * n + j] > 0) {
				matrix[i * n + j] = weight(merstw);
				matrix[j * n + i] = matrix[i * n + j];
			}
		}
	}

}



Graph::~Graph() {
	resize(0);
}
