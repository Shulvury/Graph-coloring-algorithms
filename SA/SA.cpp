#include "SA.h"
#include <random>
#include <iostream>
#include <cmath>
#include "initializeColoring.h"

using namespace std;

extern unsigned long long maxChecks;
extern unsigned long long numConfChecks;

extern double alpha; // cooling rate for temperature
extern double t_min;
extern clock_t clockStart;

uniform_real_distribution<double> unif(0.0, 1.0);
uniform_int_distribution<int> coinflip(0, 1);
random_device rd;
mt19937 mt(rd());

//this function is for modified penalty annealing
static bool TryGreedy(Graph& g, vector<int>& col, vector<vector<int>>& ColClass, vector<int>& E, int k, vector<vector<int>>& adjList,
						int& energy, int& cost) {

	int pos = -1;
	for (int i = 0; i < k; i++) {
		numConfChecks++;
		if (E[i] > 0) {
			pos = i;
			break;
		}
	}
	int node1, node2;
	for (auto v : ColClass[pos]) {
		for (auto u : ColClass[pos]) {
			numConfChecks++;
			if (u == v) continue;
			if (g[u][v] > 0) {
				node1 = v;
				node2 = u;
				goto next_stage;	//three lines below, out of this loop
			}
		}
	}
next_stage:
	//try to find color for node 1
	bool asgn = true;
	for (int c = 0; c < k; c++) {
		if (c == pos) continue;
		asgn = true;
		for (auto v : ColClass[c]) {
			if (g[v][node1] > 0) {
				asgn = false;
				break;
			}
		}
		if (asgn == true) {
			col[node1] = c;
			return true;
		}
	}
	//now do same for node 2, if node 1 failed
	for (int c = 0; c < k; c++) {
		if (c == pos) continue;
		asgn = true;
		for (auto v : ColClass[c]) {
			if (g[v][node2] > 0) {
				asgn = false;
				break;
			}
		}
		if (asgn == true) {
			col[node2] = c;
			return true;
		}
	}

	//if (RandomReassign == false) return false; //if option is disabled, we return. else, we try to randomly assign new color

	//if we are here it means we couldnt find any color with greedy approach

	//lets do random assignment then
	//before, node1 < node2. lets randomize it 
	node1 = coinflip(mt) ? node1 : node2;

	uniform_int_distribution<int> rc(0, k - 1);
	int newcol = rc(mt);
	numConfChecks++;
	while (newcol == pos) newcol = rc(mt);

	//find position of this node inside the color class
	int vert_position = -1;
	for (int v = 0; v < ColClass[pos].size(); v++) {
		if (ColClass[pos][v] == node1) vert_position = v;
	}

	int delta = SimanDeltaEnergy(g, col, ColClass, node1, newcol, E);

	for (auto v : adjList[node1]) {
		if (col[v] == pos) {
			E[pos]--;
			energy--;
		}
		if (col[v] == newcol) {
			E[newcol]++;
			energy++;
		}
	}
	ColClass[pos].erase(ColClass[pos].begin() + vert_position);
	ColClass[newcol].push_back(node1);
	col[node1] = newcol;
	numConfChecks += 2 * adjList[node1].size();

	cost = cost + delta;

	return false;
}

//this is same function but for normal annealing
static bool TryGreedy(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int k, vector<vector<int>>& M, int& E) {

	//identify problematic nodes
	int node1, node2;

	for (int i = 0; i < g.n; i++) {
		for (auto v : adjList[i]) {
			numConfChecks++;
			if (col[i] == col[v]) {
				node1 = i;
				node2 = v;
				goto next_stage;
			}
		}
	}
next_stage:
	
	bool asgn = true;
	for (int c = 0; c < k; c++) {
		numConfChecks++;
		asgn = true;
		if (c == col[node1]) continue;
		for (int v = 0; v < g.n; v++) {
			if (col[v] == c) {
				numConfChecks++;
				if (g[v][node1] > 0) {
					asgn = false;
					break;
				}
			}
		}
		numConfChecks++;
		if (asgn) {
			E = E + M[node1][c] - M[node1][col[node1]];
			for (auto i : adjList[node1]) {
				M[i][col[node1]]--;
				M[i][c]++;
			}
			col[node1] = c;
			return true;
		}
	}
	for (int c = 0; c < k; c++) {
		numConfChecks++;
		asgn = true;
		if (c == col[node2]) continue;
		for (int v = 0; v < g.n; v++) {
			if (col[v] == c) {
				numConfChecks++;
				if (g[v][node2] > 0) {
					asgn = false;
					break;
				}
			}
		}
		numConfChecks++;
		if (asgn) {
			E = E + M[node2][c] - M[node2][col[node2]];
			for (auto i : adjList[node2]) {
				M[i][col[node2]]--;
				M[i][c]++;
			}
			col[node2] = c;
			return true;
		}
	}

	//if we are here it means we were unsuccessfull and no legal recoloring is possible at all.
	// temperature is quite probably low enough to guarantee we wont be able to do "worse" move ever again
	// this means we now need to execute random move, and hopefully reach legal coloring later

	//before, node1 < node2. lets randomize it 
	node1 = coinflip(mt) ? node1 : node2;

	uniform_int_distribution<int> rc(0, k-1);
	int newcol = rc(mt);
	numConfChecks++;
	while (newcol == col[node1]) newcol = rc(mt);

	//now simply color it!
	E = E + M[node1][newcol] - M[node1][col[node1]];
	for (auto i : adjList[node1]) {
		M[i][col[node1]]--;
		M[i][newcol]++;
	}
	col[node1] = newcol;


	return false;
}

// this is same function for original simulated annealing
static bool TryGreedy(Graph& g, vector<int>& col, vector<vector<int>>& ColClass, vector<int>& E, int k, vector<vector<int>>& adjList,
	int& energy, int& cost, bool RandomReassign) {

	int pos = -1;
	for (int i = 0; i < k; i++) {
		numConfChecks++;
		if (E[i] > 0) {
			pos = i;
			break;
		}
	}
	int node1 = -1, node2 = -1;
	int nodepos1, nodepos2;
	bool found_vertices = false;
	for (int i = 0; i < ColClass[pos].size(); i++) {
		int v = ColClass[pos][i];
		if (found_vertices == true) break;
		for (int j = 0; j < ColClass[pos].size(); j++) {
			int u = ColClass[pos][j];
			numConfChecks++;
			if (u == v) continue;
			if (g[u][v] > 0) {
				node1 = v;
				nodepos1 = i;
				node2 = u;
				nodepos2 = j;
				found_vertices = true;
			}
		}
	}

	if (node1 == -1) return false;

	//try to find color for node 1
	bool asgn = true;
	for (int c = 0; c < k; c++) {
		if (c == pos) continue;
		asgn = true;
		for (auto v : ColClass[c]) {
			if (g[v][node1] > 0) {
				asgn = false;
				break;
			}
		}
		if (asgn == true) {
			int delta = SimanDeltaEnergy(g, col, ColClass, node1, c, E);
			ColClass[pos].erase(ColClass[pos].begin() + nodepos1);
			ColClass[c].push_back(node1);
			col[node1] = c;

			for (auto v : adjList[node1]) {
				if (col[v] == pos) {
					E[pos]--;
					energy--;
				}
				if (col[v] == c) {
					E[c]++;
					energy++;
				}
			}
			numConfChecks += 2 * adjList[node1].size();

			cost = cost + delta;

			return true;
		}
	}
	//now do same for node 2, if node 1 failed
	for (int c = 0; c < k; c++) {
		if (c == pos) continue;
		asgn = true;
		for (auto v : ColClass[c]) {
			if (g[v][node2] > 0) {
				asgn = false;
				break;
			}
		}
		if (asgn == true) {
			int delta = SimanDeltaEnergy(g, col, ColClass, node2, c, E);
			ColClass[pos].erase(ColClass[pos].begin() + nodepos2);
			ColClass[c].push_back(node2);
			col[node2] = c;

			for (auto v : adjList[node2]) {
				if (col[v] == pos) {
					E[pos]--;
					energy--;
				}
				if (col[v] == c) {
					E[c]++;
					energy++;
				}
			}
			numConfChecks += 2 * adjList[node2].size();

			cost = cost + delta;

			return true;
		}
	}

	return false;
}

void makeAdjList(vector< vector<int> >& adjList, Graph& g)
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

void DSaturKcol(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, int k) {

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
	candSol.push_back(vector<int>());
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
		assignAColourDSatur(foundColour, candSol, permutation, nodePos, satDeg, g, colNode, adjList);
		if (!foundColour) {
			//If we are here we have to make a new colour as we have tried all the other ones and none are suitable

			//Addendum: if we exceed needed colors, we skip and store it for later use

			if (candSol.size() >= k) {
				cout << " ";
			}
			else { // execute normal routine
				candSol.push_back(vector<int>());
				candSol.back().push_back(permutation[nodePos]);
				colNode[permutation[nodePos]] = candSol.size() - 1;
				//Remember to update the saturation degree array
				for (i = 0; i < permutation.size(); i++) {
					numConfChecks++;
					if (g[permutation[nodePos]][permutation[i]]) {
						satDeg[i]++;
					}
				}
			}
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


int SimanEnergy(Graph& g, vector<int>& col, vector<vector<int>>& ColClass, int k, vector<int>& E) {

	int cost = 0;
	for (int i = 0; i < k; i++) {
		int E = 0;
		numConfChecks += ColClass[i].size() * (ColClass[i].size() - 1);
		for (auto v : ColClass[i]) {
			for (auto u : ColClass[i]) {
				if (v == u) continue;
				if (g[v][u] > 0) E++;
			}
		}
		cost += ColClass[i].size() * (2 * E - ColClass[i].size());			// -C^2 + 2CE = C*(2E - C)
	}
	return cost;
}

int SimanDeltaEnergy(Graph& g, vector<int>& col, vector<vector<int>>& ColClass, int node, int newcol, vector<int>& E) {

	int x = ColClass[col[node]].size();
	int y = ColClass[newcol].size();

	int delta = 2 * x - 2 * y - 2;
	int E_old = 0, E_new = 0;

	for (auto v : ColClass[col[node]]) {
		if (v == node) continue;
		if (g[v][node] > 0) E_old++;
	}
	delta = delta - 2 * x * E_old - 2 * E_old - 2 * E[col[node]];

	for (auto v : ColClass[newcol]) {
		if (g[v][node] > 0) E_new++;
	}
	delta = delta + 2 * y * E_new + 2 * E_new + 2 * E[newcol];

	/*int delta2 = x * x + y * y - (x - 1) * (x - 1) - (y + 1) * (y + 1)
				- 2 * x * E[col[node]] - 2 * y * E[newcol] + 2 * (x - 1) * (E[col[node]] - E_old) + 2 * (y + 1) * (E[newcol] + E_new);

	if (delta2 != delta) {
		cout << "Achtung wrong calc in energy!!\n";
		cout << "Achtung wrong calc in energy!!\n";
		cout << "Achtung wrong calc in energy!!\n";
		cout << "Achtung wrong calc in energy!!\n";
	}*/

	return delta;
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
			for (auto v : adjList[permutation[nodePos]]) {
				if (!foundColor) break;
				for (auto t : T) {
					if (colNode[v] - i == t && colNode[v] > -1) {
						foundColor = false;
					}
					if (colNode[v] - i == -t && colNode[v] > -1) {
						foundColor = false;
					}
				}
			}
			if (foundColor) {
				colNode[permutation[nodePos]] = i;
				break;
			}
		}
		for (i = 0; i < permutation.size(); i++) {
			numConfChecks++;
			if (g[permutation[nodePos]][permutation[i]]) {
				satDeg[i]++;
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

		for (int c = 0; c < k; c++) { //for all colors
			foundColour = true;
			for (auto v : adjList[permutation[nodePos]]) { //for all adjacency list
				if (!foundColour) break;
				for (int t = 0; t < T.size(); t++) {
					if (!foundColour) break;
					if (c - colNode[v] == T[t] && colNode[v] > -1) foundColour = false; //check if assigning C to V makes any illegal difference
					if (c - colNode[v] == -T[t] && colNode[v] > -1) foundColour = false;
				}
			}
			if (foundColour) {
				candSol[c].push_back(permutation[nodePos]);
				colNode[permutation[nodePos]] = c;
				break;
			}
		}
		//update saturation array
		for (i = 0; i < permutation.size(); i++) {
			numConfChecks++;
			if (g[permutation[nodePos]][permutation[i]]) {
				satDeg[i]++;
			}
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

void SimulatedAnnealing(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList, vector< vector<int> >& M, 
						vector<vector<int>>& ColClass, double& temp, int& E, int& loop) {

	int counter = 0;
	uniform_int_distribution<int> uint_color(0, k - 1);

	int good_trig = 0;
	int bad_trig = 0;
	int changes = 0;

	int cutoff = 0.3 * g.n;
	int minpercent = 0.01 * g.n;

	int freezing = 0;

	clock_t enter_point = clock();

	int zero_changes = 0;

	/*vector<vector<int>> vert_in_conf;
	for (int i = 0; i < g.n; i++) {
		vert_in_conf.push_back({});
	}
	for (auto c : ColClass) {
		for (auto v : c) {
			for (auto u : c) {
				numConfChecks++;
				if (v == u) continue;
				if (g[v][u] > 0) {
					vert_in_conf[v].push_back(u);
					vert_in_conf[u].push_back(v);
				}
			}
		}
	}*/
	

	while (E != 0 && numConfChecks < maxChecks) {

		numConfChecks += 3;
		
		if (E == 1) {
			if (TryGreedy(g, col, adjList, k, M, E) == true) {
				E = 0;
				for (int i = 0; i < g.n; i++) {
					E += M[i][col[i]];
				}
				if (0 != E) cout << "ACHTUUUUUUUUUUUUUUUNG!\n";
				return;
			}
		}

		//select vertices that are in conflict
		vector<int> conf_vertex;
		for (int i = 0; i < g.n; i++) {
			if (M[i][col[i]] > 0) {
				conf_vertex.push_back(i);
			}
		}
		numConfChecks += g.n;
		int n = conf_vertex.size();
		uniform_int_distribution<int> uint_node(0, n - 1);

		int node = conf_vertex[uint_node(mt)]; // random vertex  
		int c = uint_color(mt);		 // new random color
		while (c == col[node]) {
			c = uint_color(mt); // make sure color is different from current one
			numConfChecks++;
		}

		//now evaluate how good move is
		int newval = M[node][c] - M[node][col[node]];
		// if its at least as good its always accepted
		if (newval <= 0) {

			if (newval == 0) zero_changes++;

			good_trig++; 
			changes++;

			//update cost matrix
			for (auto i : adjList[node]) {
				M[i][col[node]]--;
				M[i][c]++;
			}
			
			//assign new color and set new energy value
			col[node] = c;
			E = E + newval;
			numConfChecks++;
		}
		//lets see if it is accepted by probability
		else {
			double r = unif(mt);
			if (r < exp(-newval / temp)) {
				bad_trig++;
				changes++;
				//do the same
				for (auto i : adjList[node]) {
					M[i][col[node]]--;
					M[i][c]++;
				}
				col[node] = c;
				E = E + newval;
				numConfChecks++;
			}
		}
		// if we had enough runs at current temperature temp, reduce temp and start counting anew
		counter++;
		numConfChecks++;

		numConfChecks += 3;
		if (counter > loop || changes > cutoff) {
			if (changes > minpercent) {
				freezing = 0;
			}
			else freezing++;
			if (zero_changes - changes == 0) {
				bool temp_variable = TryGreedy(g, col, adjList, k, M, E);	//if we are stuck, try to randomize a bit
			}
			counter = 0;
			temp = alpha * temp;
			changes = 0;
			zero_changes = 0;
			if ((clock() - enter_point) / CLOCKS_PER_SEC > 3) {
				enter_point = clock();
				cout << "Current values: energy(edges) " << E << " || temperature: " << temp << "\n";
			}
		}
		if (freezing == 10) {
			cout << "---Reached freezing temp at " << temp << endl;
			return;
		}

	}

	cout << "During this run we observed " << good_trig << " posivitive changes, and " << bad_trig << " number of worse changes.\n\n";
}


void SpanSimulatedAnnealing(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList, double& temp, int& E, int& loop, vector<int>& T) {

	int counter = 0;
	uniform_int_distribution<int> uint_color(0, k - 1);

	while (E != 0 && numConfChecks < maxChecks && temp > t_min) {

		//select vertices that are in conflict
		vector<int> temporary_vector(g.n, 0);
		vector<int> conf_vertex;
		for (int i = 0; i < g.n; i++) {
			int alsize = adjList[i].size();
			for (auto v : adjList[i]) {
				for (int t = 0; t < T.size(); t++) {
					numConfChecks++;
					if (col[i] - col[v] == T[t] || col[i] - col[v] == -T[t]) {
						temporary_vector[i]++;
						temporary_vector[v]++;
					}
				}
			}
		}

		for (int i = 0; i < g.n; i++) {
			if (temporary_vector[i] > 0) conf_vertex.push_back(i);
		}

		//cout << "We have " << conf_vertex.size() << " conflicting vertices\n";

		int Tsize = T.size();
		//for (int i = 0; i < g.n; i++) {
		//	numConfChecks += adjList[i].size() * Tsize / 2;
		//}

		int n = conf_vertex.size();
		uniform_int_distribution<int> uint_node(0, n - 1);

		int node = conf_vertex[uint_node(mt)]; // random vertex  
		int c = uint_color(mt);		 // new random color
		while (c == col[node]) {
			c = uint_color(mt); // make sure color is different from current one
		}
		numConfChecks += 2;

		//now evaluate how good move is
		int newval = E;
		for (auto i : adjList[node]) {
			for (int t = 0; t < Tsize; t++) {
				numConfChecks += 2;
				if (col[i] - col[node] == T[t] || col[i] - col[node] == -T[t]) {
					newval--; // if we had conflict in prev color, then new color erases it
					break;
				}
				if (col[i] - c == T[t] || col[i] - c == -T[t]) {
					newval++;		   // if we got a conflict in new color, value becomes higher 
					break;
				}
			}
		}
		// if its at least as good its always accepted
		if (newval <= E) {
			//assign new color and set new energy value
			col[node] = c;
			E = newval;
			numConfChecks += 2;
		}
		//lets see if it is accepted by probability
		else {
			double r = unif(mt);
			if (r < exp((E - newval) / temp)) {
				col[node] = c;
				E = newval;
				numConfChecks += 2;
			}
		}
		// if we had enough runs at current temperature temp, reduce temp and start counting anew
		counter++;
		if (counter > loop) {
			counter = 0;
			temp = alpha * temp;
		}
	}

}


void ConflictMatrix(Graph& g, vector<int>& col, int k, vector< vector<int> >& M, vector< vector<int> >& adjList) {

	for (int i = 0; i < g.n; i++) {
		//numConfChecks += adjList[i].size();
		for (auto v : adjList[i]) {
			M[i][col[v]]++;
		}
	}
}


void CalculateEnergy(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int& E) {

	E = 0;
	for (int i = 0; i < g.n; i++) {
		numConfChecks += adjList[i].size();
		for (auto u : adjList[i]) {
			if (col[i] == col[u]) E++;
		}
	}
	E = E/2;
}

void SpanCalculateEnergy(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int& E, vector<int>& T) {

	//we have
	// 1 - number of violated constraints
	// 2 - sum of amounts of violations
	// 3 - span
	// 4 - number of frequencies used
	// 5 - largest frequency
	// 6 - largest amount of violation

	E = 0;
	numConfChecks += g.n;

	for (int i = 0; i < g.n; i++) {
		int s = adjList[i].size();
		numConfChecks += s;
		for (auto v : adjList[i]) {
			for (int t = 0; t < T.size(); t++) {
				if (col[i] - col[v] == T[t] || col[i] - col[v] == -T[t]) {
					E++; // i and v create conflict. move to next v in adjList
					break;
				}
			}
		}
	}
	E = E / 2;
}


void SpanInitialRandomCol(Graph& g, vector<int>& col, vector<int>& T, vector< vector<int>>& adjList) {

	vector<int> permutation(g.n);

	for (int i = 0; i < g.n; i++) permutation[i] = i;
	uniform_int_distribution<int> rv(0, g.n - 1);
	for (int i = 0; i < g.n; i++) {
		swap(permutation[rv(mt)], permutation[rv(mt)]);
	}

	for (int i = 0; i < g.n; i++) col[i] = -1;

	col[permutation[0]] = 0;

	int node;
	vector<int> extra_vert;
	int ColNum = 1;

	for (int i = 1; i < g.n; i++) {
		node = permutation[i];
		// try to assign into existing class
		bool foundColor = true;
		for (int c = 0; c < ColNum; c++) {
			foundColor = true;
			// check if its viable 
			for (auto u : adjList[node]) {
				for (auto t : T) {
					if (!foundColor) break;
					if ((col[u] > -1 && col[u] - c == t) || (col[u] > -1 && c - col[u] == t)) {
						foundColor = false;
						//cout << c << ':' << node << 'x' << u;
						//cout << "|| edge:" << g[u][node] << "|| color of u: " << col[u] << endl;
					}
				}
			}
			if (foundColor) {
				col[node] = c;
				//cout << "found color for " << node << " in " << c << endl;
				break;
			}
		}
		if (foundColor) continue;
		col[node] = ColNum;
		ColNum++;
		//cout << "\n\n\t\t " << ColNum << "\n\n";
	}


}



void FAP_calcMatrix(Graph& g, vector<int>& col, vector< vector<int> >& adjList, vector<vector<int>>& M, int k) {

	//this is just sum of violations
	for (int i = 0; i < g.n; i++) {
		//numConfChecks += adjList[i].size() * k;
		for (auto v : adjList[i]) {
			for (int c = 0; c < k; c++) {
				M[i][c] += max(0, g[i][v] - abs(c - col[v])); //0 if colors are far away enough; otherwise lenght of violation
			}
		}
	}

}

int FAP_energy(Graph& g, vector<int>& col, vector< vector<int> >& adjList, vector<vector<int>>& M, int k) {

	int E = 0;
	int mu1 = 1, mu2 = 0;
	for (int i = 0; i < g.n; i++) {
		E += M[i][col[i]] * mu1;
		if (M[i][col[i]] > 0) E += mu2;
	}

	numConfChecks += g.n;

	return E / 2;
}


void FAP_Execute_Move(Graph& g, vector<int>& col, vector< vector<int> >& adjList, vector<vector<int>>& M,
	int node, int newcolor, int& E, int k) {

	//first update energy
	E = E + M[node][newcolor] - M[node][col[node]];

	//next update cost matrix
	for (auto v : adjList[node]) {
		for (int c = 0; c < k; c++) {
			M[v][c] = M[v][c]																//old value
				- max(0, g[node][v] - abs(col[node] - c)				// minus what was the confl
					+ max(0, g[node][v] - abs(newcolor - c)));				// plus what now is the confl

		}
	}
	//finally update color
	col[node] = newcolor;

	//numConfChecks += adjList[node].size() * k;
}

int FAP_EnergyDifference(Graph& g, vector<int>& col, vector<vector<int>>& adjList, int node, int newcol) {

	int E = 0;
	for (auto v : adjList[node]) {
		E += max(0, g[node][v] - abs(col[v] - newcol)) - max(0, g[node][v] - abs(col[v] - col[node]));
	}

	return E;
}

void FAP_SimulatedAnnealing(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList,
	double& temp, int& E, int& loop) {

	int counter = 0;
	uniform_int_distribution<int> uint_color(0, k - 1);

	int good_trig = 0;
	int bad_trig = 0;

	int numOfTrig = 0;
	int consequtive_frozen_temp = 0;

	clock_t enter_point = clock();

	while (E > 0 && numConfChecks < maxChecks && temp > t_min) {

		numConfChecks += 3;


		if ((clock() - enter_point) / CLOCKS_PER_SEC > 5) {
			enter_point = clock();
			cout << "currently at " << numConfChecks / 1000 << "k operations!\n";
		}

		//select vertices that are in conflict
		vector<int> conf_vertex;
		for (int i = 0; i < g.n; i++) {
			for (auto v : adjList[i]) {
				numConfChecks++;
				if (g[i][v] > abs(col[i] - col[v])) {
					conf_vertex.push_back(i);
					break;
				}
			}
		}
		/*for (int i = 0; i < g.n; i++) { // for each color
			if (M[i][col[i]] > 0) conf_vertex.push_back(i);
		}*/


		//if (counter == 0) {
		//	for (auto v : conf_vertex) {
		//		cout << v << ' ';
		//	}
		//	cout << endl << endl;
		//}


		int n = conf_vertex.size();
		if (n == 0) return;
		uniform_int_distribution<int> uint_node(0, n - 1);

		int node = conf_vertex[uint_node(mt)]; // random vertex  
		int c = uint_color(mt);		 // new random color
		while (c == col[node]) {
			c = uint_color(mt); // make sure color is different from current one
			numConfChecks++;
		}

		//now evaluate how good move is
		//int newval = E + M[node][c] - M[node][col[node]];
		int diff = FAP_EnergyDifference(g, col, adjList, node, c);

		// if its at least as good its always accepted
		if (diff <= 0) {

			good_trig++;
			numOfTrig++;
			numConfChecks++;
			E += diff;
			col[node] = c;
			//execute move and update M, col, E
			//FAP_Execute_Move(g, col, adjList, M, node, c, E, k);
		}
		//lets see if it is accepted by probability
		else {
			numConfChecks += 2;
			double r = unif(mt);
			if (r < exp((-diff) / temp)) {
				bad_trig++;
				numOfTrig++;
				E += diff;
				col[node] = c;
				//do the same
				//execute move and update M, col, E
				//FAP_Execute_Move(g, col, adjList, M, node, c, E, k);
			}
		}
		// if we had enough runs at current temperature temp, reduce temp and start counting anew
		counter++;
		numConfChecks++;
		if (counter > loop) {
			counter = 0;
			temp = alpha * temp;
			if (numOfTrig > 0) {
				numOfTrig = 0;
				consequtive_frozen_temp = 0;
			}
			else {
				numOfTrig = 0;
				consequtive_frozen_temp++;
			}

			if (consequtive_frozen_temp == 10) {
				cout << "Frozen temperatuer for the 10th time. Terminating search.\n";
				return;
			}
		}


	}


	cout << "During this run we observed " << good_trig << " posivitive changes, and " << bad_trig << " number of worse changes.\n\n";

}

void FAP_K_col_Dsatur(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList, int k) {

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

		for (int c = 0; c < k; c++) { //for all colors
			foundColour = true;
			for (auto v : adjList[permutation[nodePos]]) { //for all adjacency list
				if (!foundColour) break;
				if (g[permutation[nodePos]][v] > abs(colNode[v] - c)) foundColour = false;
				// eligable colors satisfy | col_1 - col_2 | >= g[1][2]
			}
			if (foundColour) {
				candSol[c].push_back(permutation[nodePos]);
				colNode[permutation[nodePos]] = c;
				break;
			}
		}
		//update saturation array
		for (i = 0; i < permutation.size(); i++) {
			numConfChecks++;
			if (g[permutation[nodePos]][permutation[i]]) {
				satDeg[i]++;
			}
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

void FAP_Intial_DSatur(vector<int>& colNode, Graph& g, vector< vector<int> >& adjList) {

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

	int dctr = 0;

	while (!permutation.empty()) {

		//cout << "start " << dctr << '|';
		//dctr++;

		//choose the node to colour next (the rightmost node that has maximal satDegree)
		maxSat = INT_MIN;
		for (i = 0; i < satDeg.size(); i++) {
			if (satDeg[i] >= maxSat) {
				maxSat = satDeg[i];
				nodePos = i;
			}
		}
		//cout << "chosen node " << nodePos << '|';
		//now choose which colour to assign to the node
		int clr = 1e5;
		bool foundColor = true;
		for (int i = 0; i < clr; i++) {
			foundColor = true;
			for (auto v : adjList[permutation[nodePos]]) {
				if (!foundColor) break;
				if (colNode[v] == INT_MIN) continue;
				if (g[permutation[nodePos]][v] > abs(colNode[v] - i)) {
					foundColor = false;
					//cout << "color " << i << " is no good, cause of " << v << " edge is " << g[permutation[nodePos]][v] << endl;
					//cout << v << " violates color " << i << " for node " << nodePos << "|| their edge is:" << g[nodePos][v] << endl;
				}
				// eligable colors satisfy | col_1 - col_2 | >= g[1][2]
				numConfChecks += 3;
			}
			if (foundColor) {
				colNode[permutation[nodePos]] = i;
				break;
			}
		}
		//cout << "vertex " << nodePos << " got color " << colNode[permutation[nodePos]] << endl;
		for (i = 0; i < permutation.size(); i++) {
			numConfChecks++;
			if (g[permutation[nodePos]][permutation[i]]) {
				satDeg[i]++;
			}
		}

		//Finally, we remove the node from the permutation
		permutation.erase(permutation.begin() + nodePos);
		satDeg.erase(satDeg.begin() + nodePos);
	}

}


int FAP_InitialEnergy(Graph& g, vector<int>& col, vector<vector<int>>& adjList) {

	int E = 0;
	for (int i = 0; i < g.n; i++) {
		for (auto v : adjList[i]) {
			E += max(0, g[i][v] - abs(col[v] - col[i]));
		}
	}

	return E / 2;
}


double FAP_Seek_T0(Graph& g, vector<int>& col, int k, vector< vector<int> >& adjList) {

	int counter = 0;
	uniform_int_distribution<int> uint_color(0, k - 1);

	double temp = 1.0;

	int trig = 0;

	//select vertices that are in conflict
	vector<int> conf_vertex;
	for (int i = 0; i < g.n; i++) {
		for (auto v : adjList[i]) {
			numConfChecks++;
			if (g[i][v] > abs(col[i] - col[v])) {
				conf_vertex.push_back(i);
				break;
			}
		}
	}
	int n = conf_vertex.size();
	if (n == 0) {
		cout << "No conflicts found, returning!\n";
		return 1.0;
	}
	uniform_int_distribution<int> uint_node(0, n - 1);

	double acceptance = 0;

	while (acceptance < 0.9) {
		numConfChecks++;
		counter++;

		for (int i = 0; i < 100; i++) {
			int node = conf_vertex[uint_node(mt)]; // random vertex  
			int c = uint_color(mt);		 // new random color
			while (c == col[node]) {
				c = uint_color(mt); // make sure color is different from current one
				numConfChecks++;
			}

			int delta = FAP_EnergyDifference(g, col, adjList, node, c);
			if (delta <= 0) trig++;
			else if (unif(mt) < exp(-delta) / temp) trig++;
		}
		numConfChecks += 200;

		acceptance = trig / 100.0;
		if (acceptance < 0.9) temp *= 2;

	}

	return temp;
}


int PenaltySimulatedAnnealing(Graph& g, vector<int>& col, vector<vector<int>>& adjList, vector<vector<int>> ColClass,
	double& temp, int loop, int& k) {

	vector<int> E(k, 0);
	int energy = 0;
	for (int i = 0; i < g.n; i++) {
		for (auto v : adjList[i]) {
			if (col[i] == col[v]) E[col[v]]++;
		}
	}
	for (int i = 0; i < k; i++) {
		E[i] = E[i] / 2;
		energy += E[i];
	}

	if (energy == 0) {
		cout << "Found " << k << " coloring with DSatur\n\n";
		return 0;
	}

	clock_t duration = clock();

	int cost = SimanEnergy(g, col, ColClass, k, E);
	cout << "\tStarting energy is " << cost << " || Bad Edges: " << energy << "\n\n";
	int freezing = 0;

	int best_col = k;
	uniform_int_distribution<int> randCol(0, k - 1);

	int counter = 0;
	int change = 0;

	int energy_change = 0;
	//int previous_energy = energy;

	while (numConfChecks < maxChecks && energy > 0) {

		if (energy == 1) {
			if (TryGreedy(g, col, ColClass, E, k, adjList, energy, cost)) return 0;
		}
		//else if (energy == 2 && temp < 2e-322) bool tempvar = TryGreedy(g, col, ColClass, E, k, adjList, energy, cost);
		
		counter++;
		numConfChecks += 3;

		if ((clock() - duration) / CLOCKS_PER_SEC > 3) {
			duration = clock();
			cout << "Current values: Cost(function) " << cost << " || Energy(edges) " << energy << " || Colors " << E.size() << "\n\n";
		}

		int oldcol = randCol(mt);
		int newcol = randCol(mt);

		numConfChecks += 3;
		while (ColClass[oldcol].empty()) {
			oldcol = randCol(mt);
		}
		while (newcol == oldcol) {
			newcol = randCol(mt);
		}

		uniform_int_distribution<int> randVert(0, ColClass[oldcol].size() - 1);
		int pos = randVert(mt);
		int node = ColClass[oldcol][pos];

		int delta = SimanDeltaEnergy(g, col, ColClass, node, newcol, E);

		numConfChecks += 2;
		if (delta <= 0) {

			ColClass[oldcol].erase(ColClass[oldcol].begin() + pos);
			ColClass[newcol].push_back(node);
			col[node] = newcol;

			for (auto v : adjList[node]) {
				if (col[v] == oldcol) {
					E[oldcol]--;
					energy--;
				}
				if (col[v] == newcol) {
					E[newcol]++;
					energy++;
				}
			}
			numConfChecks += 2 * adjList[node].size();

			cost = cost + delta;
			change++;
			//if (energy != previous_energy) energy_change++;

		}
		else if (unif(mt) < exp(-delta / temp)) {

			ColClass[oldcol].erase(ColClass[oldcol].begin() + pos);
			ColClass[newcol].push_back(node);
			col[node] = newcol;

			for (auto v : adjList[node]) {
				if (col[v] == oldcol) {
					E[oldcol]--;
					energy--;
				}
				if (col[v] == newcol) {
					E[newcol]++;
					energy++;
				}
			}
			numConfChecks += 2 * adjList[node].size();

			cost = cost + delta;
			change++;
			//if (energy != previous_energy) energy_change++;
		}
		//previous_energy = energy;

		numConfChecks += 4;
		if (counter >= loop || change < 0.1 * loop) {		// 30% cutoff
			//if (temp < 2e-321 && energy_change == 0) {
			//	freezing++;
			//	//cout << "Freezing!!\n";
			//}
			//else freezing = 0;
			//if (freezing == 10000) {
			//	cout << "Unfreeze!\n";
			//	bool tempvar = TryGreedy(g, col, ColClass, E, k, adjList, energy, cost);
			//	freezing = 0;
			//}
			//energy_change = 0;
			counter = 0;
			change = 0;
			temp = temp * alpha;
			//if (freezing == 1000 && temp == t_min) {
			//	freezing = 0;
			//	bool temp_variable = TryGreedy(g, col, ColClass, E, k, adjList, energy, cost);
			//}
			//else energy_on_prev_temp = energy;
		}

		if (energy == 0 && ColClass[oldcol].empty()) {
			cout << "\t\tFound new best coloring with " << best_col << " colors\n";
			return energy;
		}

	}

	//if (freezing == 10) {
	//	cout << "Freezing " << freezing << " at temperature " << temp << endl;
	//}

	cout << "\tEnding one run with cost function = " << cost << " || Bad Edges: " << energy << " || temp: " << temp << "\n\n";

	return energy;
}



void randomizeLastColor(Graph& g, vector<int>& col, int k) {

	uniform_int_distribution<int> rc(0, k - 1);

	//reassign last color class into other classes, randmoly
	for (int i = 0; i < g.n; i++) {
		if (col[i] == k) col[i] = rc(mt);
	}
	numConfChecks += g.n;

}


int UnmodifiedPenaltySimulatedAnnealing(Graph& g, vector<int>& col, vector<vector<int>>& adjList, vector<vector<int>> ColClass,
	double& temp, int loop, int& k) {

	vector<int> E(k, 0);
	int energy = 0;
	for (int i = 0; i < g.n; i++) {
		numConfChecks += adjList[i].size();
		for (auto v : adjList[i]) {
			if (col[i] == col[v]) E[col[v]]++;
		}
	}
	for (int i = 0; i < k; i++) {
		E[i] = E[i] / 2;
		energy += E[i];
	}

	if (energy == 0) {
		cout << "Found " << k << " coloring with DSatur\n\n";
		return 0;
	}

	clock_t duration = clock();

	int cost = SimanEnergy(g, col, ColClass, k, E);
	cout << "\tStarting energy is " << cost << " || Bad Edges: " << energy << "\n\n";
	int freezing = 0;

	int lowest_col = k;
	int highest_col = k;
	int current_colors = k;
	

	int counter = 0;
	int change = 0;
	//int energy_change = 0;
	//int previous_energy = energy;

	while (numConfChecks < maxChecks) {

		counter++;

		numConfChecks += 3;
		if (energy == 1) {
			bool temp_var = TryGreedy(g, col, ColClass, E, current_colors, adjList, energy, cost, false);
		}
		if (current_colors < lowest_col && energy == 0) {
			//cout << "start ";
			lowest_col = current_colors;

			//extra caution just in case
			for (int i = 0; i < g.n; i++) {
				for (auto v : adjList[i]) {
					if (col[i] == col[v]) cout << "FAILURE!!!\n";
				}
			}

			cout << "Successfully found new best coloring with " << lowest_col << endl;
			int time_expended = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
			cout << "Time taken: " << time_expended << " || temp: " << temp << "\n";
			//t_min = temp;
			cout << "Currently at " << numConfChecks / 1000 << "k checks\n\n";
		}

		if ((clock() - duration) / CLOCKS_PER_SEC > 3) {
			duration = clock();
			cout << "Current values: Cost(function) " << cost << " || Energy(edges) " << energy 
				 << " || Colors " << E.size() << " || Checks " << numConfChecks/1000000 << "m\n\n";
		}

		uniform_int_distribution<int> oldRandCol(0, current_colors - 1);
		uniform_int_distribution<int> newRandCol(0, current_colors);
		int oldcol = oldRandCol(mt);
		int newcol = newRandCol(mt);

		numConfChecks += 3;
		while (ColClass[oldcol].empty()) {			// if we somehow sampled empty col, it will cause error in program. need to reroll
			oldcol = oldRandCol(mt);
		}
		while (newcol == oldcol) {					// need new DISTINCT color 
			newcol = newRandCol(mt);
		}
		if (newcol == current_colors) { //it means new empty color was created, so we need to create more sctuctures
			E.push_back(0);				// no vertices == no conflicts
			ColClass.push_back({});		// new col class
			current_colors++;			// obviously

			if (current_colors > highest_col) {		//just for statistics
				highest_col = current_colors;
				//cout << "\t Reached new highest color usage at " << highest_col << endl;
			}
		}

		uniform_int_distribution<int> randVert(0, ColClass[oldcol].size() - 1);
		int pos = randVert(mt);
		int node = ColClass[oldcol][pos];

		int delta = SimanDeltaEnergy(g, col, ColClass, node, newcol, E);

		numConfChecks += 2;
		if (delta <= 0) {

			ColClass[oldcol].erase(ColClass[oldcol].begin() + pos);
			ColClass[newcol].push_back(node);
			col[node] = newcol;

			for (auto v : adjList[node]) {
				if (col[v] == oldcol) {
					E[oldcol]--;
					energy--;
				}
				if (col[v] == newcol) {
					E[newcol]++;
					energy++;
				}
			}
			numConfChecks += 2 * adjList[node].size();

			cost = cost + delta;
			change++;
			//if (energy != previous_energy) energy_change++;

			//if old color is now empty, we have to delete it and update all structures
			numConfChecks++;
			if (ColClass[oldcol].empty()) {
				//cout << "Begin of erase\n";
				E.erase(E.begin() + oldcol);
				ColClass.erase(ColClass.begin() + oldcol);
				current_colors--;
				numConfChecks += g.n;
				for (int i = 0; i < g.n; i++) {		// there is no vert with col[v] = oldcol - 1; so all vertices above that go 1 color below
					if (col[i] >= oldcol) col[i]--;
				}
				//cout << "end of erase\n";
				if (E.size() != ColClass.size() || E.size() != current_colors) cout << "We have a problem!!\n";
			}
		}
		else if (unif(mt) < exp(-delta / temp)) {

			ColClass[oldcol].erase(ColClass[oldcol].begin() + pos);
			ColClass[newcol].push_back(node);
			col[node] = newcol;

			for (auto v : adjList[node]) {
				if (col[v] == oldcol) {
					E[oldcol]--;
					energy--;
				}
				if (col[v] == newcol) {
					E[newcol]++;
					energy++;
				}
			}
			numConfChecks += 2 * adjList[node].size();

			cost = cost + delta;
			change++;
			//if (energy != previous_energy) energy_change++;
			
			//if old color is now empty, we have to delete it and update all structures
			numConfChecks++;
			if (ColClass[oldcol].empty()) {
				//cout << "Begin of erase\n";
				E.erase(E.begin() + oldcol);
				ColClass.erase(ColClass.begin() + oldcol);
				current_colors--;
				numConfChecks += g.n;
				for (int i = 0; i < g.n; i++) {		// there is no vert with col[v] = oldcol - 1; so all vertices above that go 1 color below
					if (col[i] >= oldcol) col[i]--;
				}
				//cout << "end of erase\n";
				if (E.size() != ColClass.size() || E.size() != current_colors) cout << "We have a problem!!\n";
			}
		}
		//previous_energy = energy;

		numConfChecks += 3;
		if (counter >= loop) {		// cuttoff disabled ->         || change < 0.3 * loop
			if (change < 0.02 * loop) freezing++;
			else freezing = 0;
			counter = 0;
			change = 0;
			temp = temp * alpha;

			for (int C = 0; C < ColClass.size(); C++) {
				if (ColClass[C].empty()) {
					//cout << "Begin of erase\n";
					E.erase(E.begin() + C);
					ColClass.erase(ColClass.begin() + C);
					current_colors--;
					numConfChecks += g.n;
					for (int i = 0; i < g.n; i++) {		// there is no vert with col[v] = oldcol - 1; so all vertices above that go 1 color below
						if (col[i] >= C) col[i]--;
					}
					//cout << "end of erase\n";
					if (E.size() != ColClass.size() || E.size() != current_colors) cout << "We have a problem!!\n";
				}
			}


		}
	}

	if (freezing == 10) {
		cout << "Frozen! \n";
		cout << "Current values: Cost(function) " << cost << " || Energy(edges) " << energy << " || Colors " << E.size() << "\n\n";
		return energy;
	}
	
	cout << "Run of allocated resources, terminating.\n";
	cout << "Current values: Cost(function) " << cost << " || Energy(edges) " << energy << " || Colors " << E.size() << "\n\n";
	return energy;
}



void RandomCol(Graph& g, vector<int>& col, int k) {

	uniform_int_distribution<int> randomCol(0, k-1);
	for (int i = 0; i < g.n; i++) {
		col[i] = randomCol(mt);
	}

}