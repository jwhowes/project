#include <iostream>
#include <math.h>
#include <vector>
#include <array>

using namespace std;

const int NUM_VERTICES = 10;

int adj_matrix[NUM_VERTICES][NUM_VERTICES];

vector<vector<int>> s;
int c;

void move(int v, int c1, int c2, vector<int> * partition) {
	// Moves vertex v from colour group c1 to colour group c2 in partition
	// Requires that v be a member of c1
	for (auto it = partition[c1].begin(); it != partition[c1].end(); ++it) {
		if (*it == v) {
			partition[c1].erase(it);
		}
	}
	partition[c2].push_back(v);
}

int f(vector<vector<int>> & partition) {
	int ret = 0;
	for (int i = 0; i < partition.size(); i++) {
		ret -= partition[i].size() * partition[i].size();
	}
	return ret;
}

bool valid(int v, int c, vector<vector<int>> & partition) {
	for (int i = 0; i < partition[c].size(); i++) {
		if (adj_matrix[v][partition[c][i]] == 1) {
			return false;
		}
	}
	return true;
}

int order[NUM_VERTICES];
void generate_initial_solution() {
	random_shuffle(begin(order), end(order));
	for (int v : order) {
		int c = 0;
		while (true) {
			if (c == s.size()) {
				s.push_back({ c });
				break;
			}
			if (valid(v, c, s)) {
				s[v].push_back(c);
				break;
			}
			c++;
		}
	}
}

vector<int> & kempe_chain(int c, int d, int v) {  // TODO
	vector<int> K = { v };
	int i = 0;
	while (i < K.size()) {
		for (int j = 0; j < NUM_VERTICES; j++) {
			if (find(K.begin(), K.end(), j) == K.end() && adj_matrix[K[i]][j] == 1) {  // j not in K

			}
		}
	}
}

int main(){
	for (int i = 0; i < NUM_VERTICES; i++) {
		order[i] = i;
	}
    std::cout << "Hello World!\n";
	return 0;
}