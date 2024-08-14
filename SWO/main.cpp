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
#include <algorithm>
#include <random>

//This makes sure the compiler uses _strtoui64(x, y, z) with Microsoft Compilers, otherwise strtoull(x, y, z) is used
#ifdef _MSC_VER
#define strtoull(x, y, z) _strtoui64(x, y, z)
#endif

using namespace std;

unsigned long long numConfChecks;

void logAndExitNow(string s) {
	//Writes message s to screen and log file and then exits the program
	ofstream resultsLog("resultsLog.log", ios::app);
	resultsLog << s;
	std::cout << s;
	resultsLog.close();
	exit(1);
}


void greedySWOcoloring(Graph& g, vector<int>& col, vector<int>& pos) {


	for (int i = 0; i < g.n; i++) {
		col[i] = 0;
	}

	numConfChecks += g.n;

	vector< vector<int> > colClasses(1);
	col[pos[0]] = 1;

	colClasses[0].push_back(pos[0]);

	numConfChecks += 2;

	bool col_is_fit = true;
	bool assigned = false;

	numConfChecks += 2;

	//for every vertex in sequence
	for (int i = 1; i < g.n; i++) {
		//check if it fits in existing classes
		for (int j = 0; j < colClasses.size(); j++) {

			if (assigned) break;
			col_is_fit = true;

			for (int l = 0; l < colClasses[j].size(); l++) {

				if (g[pos[i]][colClasses[j][l]] != 0) {
					col_is_fit = false;
					numConfChecks++;
				}
				if (!col_is_fit) break;
			}

			//assign in first fit class
			if (col_is_fit) {
				col[pos[i]] = j;
				colClasses[j].push_back(pos[i]);
				assigned = true;
				break;
			}
		}


		//if no color found, create new color
		if (!assigned) {
			numConfChecks += 3;
			colClasses.push_back(vector<int>());
			colClasses[colClasses.size() - 1].push_back(pos[i]);
			col[pos[i]] = colClasses.size() - 1;
		}

		assigned = false;

	}

}

void analyzeBlame_WithRandom(Graph& g, vector<int>& col, vector<int>& pos, int coltarget) {

	int numcol = *max_element(col.begin(), col.end());

	vector<int> blame(g.n);

	//assign blame to vertices with big color

	double p = 0.25; // top percentage
	double cof = g.n / 20.0; // maximum lenght of move

	int badcolors = numcol - coltarget;
	// 25% of top colors but at least 3 colors are "bad"
	//if (p * numcol < 3)   badcolors = 2;
	//else				  badcolors = int(p * numcol);

	for (int i = 0; i < g.n; i++) {
		if (col[i] > coltarget) {
			blame[i] = col[i] - coltarget;
		}
	}
	numConfChecks += g.n;

	double move_lenght;
	for (int i = 0; i < g.n; i++) {
		move_lenght = cof * blame[i];
		//if move exceeds size of vector, move to start
		if (int(move_lenght) > i) {
			numConfChecks++;
			for (int j = 0; j < i; j++) {
				swap(pos[j], pos[i]);
			}
		}
		//if not, move it up for how much is needed
		else {
			for (int j = move_lenght; j < i; j++) {
				swap(pos[j], pos[i]);
			}
		}
	}

	for (int i = 0; i < 2; i++) {
		swap( pos[rand() % g.n], pos[rand() % g.n] );
	}

}


void analyzeBlame_NoRandom(Graph& g, vector<int>& col, vector<int>& pos, int coltarget) {

	int numcol = *max_element(col.begin(), col.end());

	vector<int> blame(g.n);

	//assign blame to vertices with big color

	double p = 0.15; // top percentage
	double cof = g.n / 20.0; // maximum lenght of move

	int badcolors = numcol - coltarget;
	// 25% of top colors but at least 3 colors are "bad"
	//if (p * numcol < 3)   badcolors = 2;
	//else				  badcolors = int(p * numcol);

	for (int i = 0; i < g.n; i++) {
		if (col[i] > coltarget) {
			blame[i] = col[i] - coltarget;
		}
	}
	numConfChecks += g.n;

	double move_lenght;
	for (int i = 0; i < g.n; i++) {
		move_lenght = cof * blame[i];
		//if move exceeds size of vector, move to start
		if (int(move_lenght) > i) {
			for (int j = 0; j < i; j++) {
				swap(pos[j], pos[i]);
			}
		}
		//if not, move it up for how much is needed
		else {
			for (int j = move_lenght; j < i; j++) {
				swap(pos[j], pos[i]);
			}
		}
		numConfChecks += 2 * i;
	}

}

bool pair_comp(pair<int, int> a, pair<int, int> b) {
	numConfChecks += 1;
	return (a.second > b.second);
}

void SortBlame(Graph& g, vector<int>& col, vector<int>& pos) {

	int numcol = *max_element(col.begin(), col.end()) + 1;
	numConfChecks += g.n;

	vector<int> blame(g.n);

	//assign blame to vertices with big color

	for (int i = 0; i < g.n; i++) {
		blame[i] = col[i];
	}
	numConfChecks += g.n;

	vector <pair <int, int>> node_blame(g.n);
	for (int i = 0; i < g.n; i++) {
		node_blame[i] = make_pair(i, blame[i]);
	}

	sort(node_blame.begin(), node_blame.end(), pair_comp);
	for (int i = 0; i < g.n; i++) {
		pos[i] = node_blame[i].first;
	}

}

void BlameNeighb(Graph& g, vector<int>& col, vector<vector<int>>& adjlist, vector<int>& pos, int targetcol) {

	int numcol = *max_element(col.begin(), col.end()) + 1;
	numConfChecks += g.n;

	vector<int> blame(g.n);

	int blm = 0;
	for (int i = 0; i < g.n; i++) {
		for (int j = 0; j < adjlist[i].size(); j++) {
			if (col[j] >= targetcol) blm += 1;
		}
		blame[i] = blm;
		blm = 0;
		numConfChecks += adjlist[i].size();
	}

	
	double cof = g.n / 30.0; // maximum lenght of move
	double move_lenght;
	for (int i = 0; i < g.n; i++) {
		move_lenght = cof * blame[i];
		//if move exceeds size of vector, move to start
		if (int(move_lenght) > i) {
			for (int j = 0; j < i; j++) {
				swap(pos[j], pos[i]);
			}
		}
		//if not, move it up for how much is needed
		else {
			for (int j = move_lenght; j < i; j++) {
				swap(pos[j], pos[i]);
			}
		}
		numConfChecks += 2 * i;
	}
}


bool check_color_for_conflict(Graph& g, int& node, vector<int>& col) {
	if (col.empty()) return true;
	for (int i = 0; i < col.size(); i++) {
		numConfChecks += 1;
		//cout << col[i] << ' ';
		if (g[node][col[i]] != 0) return false;
	}
	return true;
}

void JoslinGreedy(Graph& g, vector<int>& col, vector<int>& pos, int numcol, vector<vector<int>> adjList) {

	vector< vector<int> > ColClass;
	for (int i = 0; i < numcol; i++) {
		ColClass.push_back({});
	}

	vector<int> prevcol(col);
	for (int i = 0; i < g.n; i++) col[i] = -1;
	//numConfChecks += g.n + g.n + numcol;

	bool checker = true;

	for (auto i : pos) {

		checker = true;
		//cout << "node " << i << " had prev color " << prevcol[i] << endl;

		//if prev color is less than target number
		if (prevcol[i] < numcol) {
			// check if prev color can be assigned
			checker = check_color_for_conflict(g, i, ColClass[prevcol[i]]);
			//for (int j = 0; j < ColClass[prevcol[i]].size(); j++) {
			//	if (g[i][j] != 0) checker = false;
			//	numConfChecks += 1;
			//}
			numConfChecks += 2;
			if (checker) {
				ColClass[prevcol[i]].push_back(i);
				col[i] = prevcol[i];
				//cout << "first cond \n";
				continue;
			}
		}
		
		//cout << "one" << endl;
		// now try to assign greedy to color < numcol
		/*for (int j = 0; j < numcol; j++) {
			//cout << j << ' ';
			if (j == prevcol[i]) continue; //already checked
			if(!ColClass[j].empty()) checker = check_color_for_conflict(g, i, ColClass[j]);
			//cout << j << ' ';
			//for (int l = 0; l < ColClass[j].size(); l++) {
			//	if (g[i][l] != 0) checker = false;
			//}
			if (checker) {
				col[i] = j;
				ColClass[j].push_back(i);
				//cout << "second cond \n";
				break;
			}
		}*/

		// try to assign "least constraining color"
		vector<int> good_colors;
		for (int c = 0; c < numcol; c++) {
			//see if its a feasible color
			if (ColClass[c].empty()) continue;
			else if (check_color_for_conflict(g, i, ColClass[c])) good_colors.push_back(c);
		}
		//now choose best
		if (!good_colors.empty()) {
			numConfChecks++;
			int measure_minimal = INT_MAX;
			int at_minimal = 0;
			int temp_measure = 0;

			int best_col = -1;
			int best_mes = -1;
			int best_at_min = INT_MAX;
			for (int j = 0; j < good_colors.size(); j++) {
				int c = good_colors[j];
				//assume c is color we chose for current node i
				for (auto v : adjList[i]) {
					temp_measure = 0;
					//see how many different colors node v can have with col[i] = c
					for (int temp_color = 0; temp_color < numcol; temp_color++) {
						numConfChecks++;
						if (temp_color == c) continue;
						if (ColClass[temp_color].empty()) temp_measure++;
						else if (check_color_for_conflict(g, v, ColClass[temp_color])) temp_measure++;
					}
					if (measure_minimal > temp_measure) {
						measure_minimal = temp_measure;
						at_minimal = 1;
						numConfChecks++;
					}
					else if (measure_minimal == temp_measure) {
						at_minimal++;
						numConfChecks += 2;
					}
				}
				if (measure_minimal > best_mes) {
					numConfChecks++;
					best_mes = measure_minimal;
					best_col = c;
					best_at_min = at_minimal;
				}
				else if (measure_minimal == best_mes && best_at_min > at_minimal) {
					numConfChecks += 2;
					best_mes = measure_minimal;
					best_col = c;
					best_at_min = at_minimal;
				}
			}

			col[i] = best_col;
			ColClass[best_col].push_back(i);

			continue;
		}

		//cout << "two" << endl;
		


		// try to grab existing color < numcol

		//check if there are any colors to grab
		vector< pair<int, int> > grab;
		for (int j = 0; j < numcol; j++) {

			vector<int> conflicting;
			// see what is in conflict with new vert 'i' in color 'j'
			for (auto v : ColClass[j]) {
				if (g[v][i] != 0) conflicting.push_back(v);
			}
			numConfChecks += ColClass[j].size();
			//check if they all can be reassigned

			bool reasgn = false; //keep track if vertex can be re-assigned somewhere
			for (auto v : conflicting) {
				//look through all colors
				reasgn = false;
				for (int k = 0; k < numcol; k++) {
					//skip its own col
					if (k == j) continue;
					numConfChecks++;
					if (check_color_for_conflict(g, v, ColClass[k])) {
						// this conflicting vertex found at least one way out
						reasgn = true;
						break;
					}
				}
				// if some vert is not reasgn, color cannot be grabbed
				numConfChecks++;
				if (!reasgn) break;
			}
			//save col and amount of vert to be reassigned 
			if(reasgn) grab.push_back(make_pair(j, conflicting.size()));
			numConfChecks++;
		}
		//if we found colors to grab
		//cout << "three" << endl;
		numConfChecks++;
		if (!grab.empty()) {

			//find least contraining color
			pair<int, int> t(-1, g.n);
			for (auto p : grab) {
				if (t.second > p.second) t = p;
			}
			numConfChecks += grab.size();

			//initialize confl vector again
			vector<int> conflicting;
			// see what is in conflict
			for (auto v : ColClass[t.first]) {
				if (g[v][i] != 0) conflicting.push_back(v);
			}
			numConfChecks += ColClass[t.first].size();
			//now find new colors for all of them
			for (auto v : conflicting) { // all confl vert
				for (int k = 0; k < numcol; k++) { //all colors except its own
					if (k == t.first) continue;
					if (check_color_for_conflict(g, v, ColClass[k])) {
						col[v] = k;
						ColClass[t.first].erase( find(ColClass[t.first].begin(), ColClass[t.first].end(), v) );
						ColClass[k].push_back(v);
						break;
					}
				}

			}
			//finally assign color to 'i'
			col[i] = t.first;
			ColClass[t.first].push_back(i);

			//all is done, now proceed with next vertex
			continue;
		}
		//cout << "four" << endl;
		// we are here if 'i' couldnt be assigned anywhere among first numcol colors
		//now just greedily assign it anywhere

		//see if we have another color for that

		for (int k = numcol; k < ColClass.size(); k++) {
			if (check_color_for_conflict(g, i, ColClass[k])) {
				col[i] = k;
				ColClass[k].push_back(i);
				break;
			}
		}

		//no existing color can take i
		//create new color and leave it there
		ColClass.push_back({i});
		col[i] = ColClass.size() - 1;
		
		//cout << "now node " << i << " has color " << col[i] << endl;
	}

}



void JoslinGreedy_OLD(Graph& g, vector<int>& col, vector<int>& pos, int numcol) {

	vector< vector<int> > ColClass;
	for (int i = 0; i < numcol; i++) {
		ColClass.push_back({});
	}

	vector<int> prevcol(col);
	for (int i = 0; i < g.n; i++) col[i] = 0;
	numConfChecks += g.n + g.n + numcol;

	bool checker = true;

	for (auto i : pos) {

		checker = true;
		//cout << "node " << i << " had prev color " << prevcol[i] << endl;

		//if prev color is less than target number
		if (prevcol[i] < numcol) {
			// check if prev color can be assigned
			checker = check_color_for_conflict(g, i, ColClass[prevcol[i]]);
			//for (int j = 0; j < ColClass[prevcol[i]].size(); j++) {
			//	if (g[i][j] != 0) checker = false;
			//	numConfChecks += 1;
			//}
			if (checker) {
				ColClass[prevcol[i]].push_back(i);
				col[i] = prevcol[i];
				//cout << "first cond \n";
				continue;
			}
		}
		//cout << "one" << endl;
		// now try to assign greedy to color < numcol
		for (int j = 0; j < numcol; j++) {
			//cout << j << ' ';
			if (j == prevcol[i]) continue; //already checked
			if (!ColClass[j].empty()) checker = check_color_for_conflict(g, i, ColClass[j]);
			//cout << j << ' ';
			//for (int l = 0; l < ColClass[j].size(); l++) {
			//	if (g[i][l] != 0) checker = false;
			//}
			if (checker) {
				col[i] = j;
				ColClass[j].push_back(i);
				//cout << "second cond \n";
				break;
			}
		}
		//cout << "two" << endl;
		// now try to grab existing color < numcol

		//check if there are any colors to grab
		vector< pair<int, int> > grab;
		for (int j = 0; j < numcol; j++) {

			vector<int> conflicting;
			// see what is in conflict with new vert 'i' in color 'j'
			for (auto v : ColClass[j]) {
				if (g[v][i] != 0) conflicting.push_back(v);
			}
			//check if they all can be reassigned

			bool reasgn = false; //keep track if vertex can be re-assigned somewhere
			for (auto v : conflicting) {
				//look through all colors
				reasgn = false;
				for (int k = 0; k < numcol; k++) {
					//skip its own col
					if (k == j) continue;

					if (check_color_for_conflict(g, v, ColClass[k])) {
						// this conflicting vertex found at least one way out
						reasgn = true;
						break;
					}
				}
				// if some vert is not reasgn, color cannot be grabbed
				if (!reasgn) break;
			}
			//save col and amount of vert to be reassigned 
			if (reasgn) grab.push_back(make_pair(j, conflicting.size()));
		}
		//if we found colors to grab
		//cout << "three" << endl;
		if (!grab.empty()) {

			//find least contraining color
			pair<int, int> t(-1, g.n);
			for (auto p : grab) {
				if (t.second > p.second) t = p;
			}

			//initialize confl vector again
			vector<int> conflicting;
			// see what is in conflict
			for (auto v : ColClass[t.first]) {
				if (g[v][i] != 0) conflicting.push_back(v);
			}
			//now find new colors for all of them
			for (auto v : conflicting) { // all confl vert
				for (int k = 0; k < numcol; k++) { //all colors except its own
					if (k == t.first) continue;
					if (check_color_for_conflict(g, v, ColClass[k])) {
						col[v] = k;
						ColClass[t.first].erase(find(ColClass[t.first].begin(), ColClass[t.first].end(), v));
						ColClass[k].push_back(v);
						break;
					}
				}

			}
			//finally assign color to 'i'
			col[i] = t.first;
			ColClass[t.first].push_back(i);

			//all is done, now proceed with next vertex
			continue;
		}
		//cout << "four" << endl;
		// we are here if 'i' couldnt be assigned anywhere among first numcol colors
		//now just greedily assign it anywhere

		//see if we have another color for that
		bool asgn = false; //keep track is 'i' is assigned

		for (int k = numcol; k < ColClass.size(); k++) {
			if (check_color_for_conflict(g, i, ColClass[k])) {
				col[i] = k;
				ColClass[k].push_back(i);
				asgn = true;
				break;
			}
		}
		//no existing color can take i
		//create new color and leave it there
		if (!asgn) {
			ColClass.push_back({ i });
			col[i] = ColClass.size() - 1;
		}

		//cout << "now node " << i << " has color " << col[i] << endl;

	}



}

void SpanJoslinGreedy(Graph& g, vector<int>& col, vector<int>& pos, int numcol) {

	vector< vector<int> > ColClass;
	for (int i = 0; i < numcol; i++) {
		ColClass.push_back({});
	}

	vector<int> prevcol(col);
	for (int i = 0; i < g.n; i++) col[i] = 0;
	numConfChecks += g.n + g.n + numcol;

	bool checker = true;

	for (auto i : pos) {

		checker = true;
		//cout << "node " << i << " had prev color " << prevcol[i] << endl;

		//if prev color is less than target number
		if (prevcol[i] < numcol) {
			// check if prev color can be assigned
			checker = check_color_for_conflict(g, i, ColClass[prevcol[i]]);
			//for (int j = 0; j < ColClass[prevcol[i]].size(); j++) {
			//	if (g[i][j] != 0) checker = false;
			//	numConfChecks += 1;
			//}
			if (checker) {
				ColClass[prevcol[i]].push_back(i);
				col[i] = prevcol[i];
				//cout << "first cond \n";
				continue;
			}
		}
		//cout << "one" << endl;
		// now try to assign greedy to color < numcol
		for (int j = 0; j < numcol; j++) {
			//cout << j << ' ';
			if (j == prevcol[i]) continue; //already checked
			if (!ColClass[j].empty()) checker = check_color_for_conflict(g, i, ColClass[j]);
			//cout << j << ' ';
			//for (int l = 0; l < ColClass[j].size(); l++) {
			//	if (g[i][l] != 0) checker = false;
			//}
			if (checker) {
				col[i] = j;
				ColClass[j].push_back(i);
				//cout << "second cond \n";
				break;
			}
		}
		//cout << "two" << endl;
		// now try to grab existing color < numcol

		//check if there are any colors to grab
		vector< pair<int, int> > grab;
		for (int j = 0; j < numcol; j++) {

			vector<int> conflicting;
			// see what is in conflict with new vert 'i' in color 'j'
			for (auto v : ColClass[j]) {
				if (g[v][i] != 0) conflicting.push_back(v);
			}
			//check if they all can be reassigned

			bool reasgn = false; //keep track if vertex can be re-assigned somewhere
			for (auto v : conflicting) {
				//look through all colors
				reasgn = false;
				for (int k = 0; k < numcol; k++) {
					//skip its own col
					if (k == j) continue;

					if (check_color_for_conflict(g, v, ColClass[k])) {
						// this conflicting vertex found at least one way out
						reasgn = true;
						break;
					}
				}
				// if some vert is not reasgn, color cannot be grabbed
				if (!reasgn) break;
			}
			//save col and amount of vert to be reassigned 
			if (reasgn) grab.push_back(make_pair(j, conflicting.size()));
		}
		//if we found colors to grab
		//cout << "three" << endl;
		if (!grab.empty()) {

			//find least contraining color
			pair<int, int> t(-1, g.n);
			for (auto p : grab) {
				if (t.second > p.second) t = p;
			}

			//initialize confl vector again
			vector<int> conflicting;
			// see what is in conflict
			for (auto v : ColClass[t.first]) {
				if (g[v][i] != 0) conflicting.push_back(v);
			}
			//now find new colors for all of them
			for (auto v : conflicting) { // all confl vert
				for (int k = 0; k < numcol; k++) { //all colors except its own
					if (k == t.first) continue;
					if (check_color_for_conflict(g, v, ColClass[k])) {
						col[v] = k;
						ColClass[t.first].erase(find(ColClass[t.first].begin(), ColClass[t.first].end(), v));
						ColClass[k].push_back(v);
						break;
					}
				}

			}
			//finally assign color to 'i'
			col[i] = t.first;
			ColClass[t.first].push_back(i);

			//all is done, now proceed with next vertex
			continue;
		}
		//cout << "four" << endl;
		// we are here if 'i' couldnt be assigned anywhere among first numcol colors
		//now just greedily assign it anywhere

		//see if we have another color for that
		bool asgn = false; //keep track is 'i' is assigned

		for (int k = numcol; k < ColClass.size(); k++) {
			if (check_color_for_conflict(g, i, ColClass[k])) {
				col[i] = k;
				ColClass[k].push_back(i);
				asgn = true;
				break;
			}
		}
		//no existing color can take i
		//create new color and leave it there
		if (!asgn) {
			ColClass.push_back({ i });
			col[i] = ColClass.size() - 1;
		}

		//cout << "now node " << i << " has color " << col[i] << endl;

	}


}


vector<int> calculateSequence(Graph& g, vector<int>& col) {

	vector<int> position;

	int numcol = *max_element(col.begin(), col.end()) + 1;
	numConfChecks += g.n;

	for (int j = 0; j <= numcol; j++) {
		for (int i = 0; i < g.n; i++) {
			if (col[i] == j){
				position.push_back(i);
			}
		}
	}
	numConfChecks += g.n;

	return position;
}

void Fill_adj_list(Graph& g, vector< vector<int> >& adj) {
	
	for (int i = 0; i < g.n; i++) {
		for (int j = 0; j < g.n; j++) {
			if (i == j) continue;
			if (g[i][j] != 0) {
				adj[i].push_back(j);
			}
		}
	}
}


int main(int argc, char** argv)
{
	Graph g;
	int k, frequency = 0, increment = 0, verbose = 0, randomSeed = 2, tenure = 0, algorithm = 1, cost, duration, constructiveAlg = 1, targetCols = 1;
	unsigned long long maxChecks = 1e9;

	tenure++;
	verbose++;
	algorithm = 2;
	inputDimacsGraph(g, "dsjr500.1c.col.txt");

	int random_mode = 0;
	// 0 == no random
	// 1 == with random

	int old_version = 1;
	// 0 == new version used
	// 1 == old version used


	random_device rd;
	mt19937 mt(rd());
	uniform_int_distribution<int> rS(0, 1e7);
	randomSeed = rS(mt);

	//This variable keeps count of the number of times information about the instance is looked up 
	numConfChecks = 0;

	//Seed
	srand(randomSeed);

	//ofstream fout("resultsLog.log", ios::app);

	//Make the adjacency list structure 
	//int** neighbors = new int* [g.n];
	//makeAdjList(neighbors, g);

	//The solution is held in the following array
	//int* coloring = new int[g.n];
	int* bestColouring = new int[g.n];

	vector< vector<int> > adjlist;
	for (int i = 0; i < g.n; i++) {
		adjlist.push_back({});
	}
	Fill_adj_list(g, adjlist);

	//Now start the timer
	clock_t clockStart = clock();

	vector<int> colouring(g.n);
	vector<int> pos(g.n);


	//Generate the initial value for k using greedy or dsatur algorithm
	k = generateInitialK(g, constructiveAlg, bestColouring);
	
	for (int i = 0; i < g.n; i++) {
		colouring[i] = bestColouring[i];
	}

	int numcol;

	cout << "Starting with " << k << " colors.\n";
	pos = calculateSequence(g, colouring);
	//MAIN ALGORITHM
	//k--;
	int counter = 0;
	while (numConfChecks < maxChecks) {

		counter++;

		//for (auto i : colouring) {
		//	cout << i << ' ';
		//}
		//cout << "\n\n --- \n\n";
		//for (auto i : pos) {
		//	cout << i << ' ';
		//}


		// k = best number of colors, enumerated 1 to k
		// k - 1 = biggest number inside col vector
		// k - 2 = next target in optimization
		pos = calculateSequence(g, colouring);

		if(random_mode == 0) analyzeBlame_NoRandom(g, colouring, pos, k - 2);
		else {
			analyzeBlame_WithRandom(g, colouring, pos, k - 2);
		}
		
		//SortBlame(g, colouring, pos);
		//BlameNeighb(g, colouring, adjlist, pos, k - 2);


		//greedySWOcoloring(g, colouring, pos);
		
		if(old_version == 0) JoslinGreedy(g, colouring, pos, k - 2, adjlist);
		else {
			JoslinGreedy_OLD(g, colouring, pos, k - 2);
		}
		

		//for (int i = 0; i < g.n; i++) {
		//	coloring[i] = colouring[i];
		//}

		//tabu(g, coloring, k, maxChecks, tenure, verbose, frequency, increment, neighbors);

		//Algorithm has finished at this k

		numcol = *max_element(colouring.begin(), colouring.end()) + 1;
		numConfChecks += g.n;

		if(counter % 1000 == 0) cout << "Stage " << counter << " has " << numcol << " colors\n";

		if (k > numcol) {
			for (int i = 0; i < g.n; i++) {
				bestColouring[i] = colouring[i];
			}
			k = numcol;
			for (int i = 0; i < g.n; i++) {
				for (int j = 0; j < g.n; j++) {
					if (g[i][j] && colouring[i] == colouring[j]) cout << "Achtung!\n";
				}
			}
			cout << "Found new best col using " << k << ", spent " << numConfChecks/1000 << "k operations\n";
			duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);
			cout << "Time spent: " << duration << "\n\n";
			numConfChecks += g.n;
		}

		//Decrement k (if the run time hasn't been reached, we'll carry on with this new value)
		//k--;
	}

	duration = int(((double)(clock() - clockStart) / CLOCKS_PER_SEC) * 1000);


	//Delete the arrays and end	
	//delete[] coloring;
	delete[] bestColouring;
	//for (int i = 0; i < g.n; i++) delete[] neighbors[i];
	//delete[] neighbors;

	
	//fout.close();



	return 0;
}



