#define _SECURE_SCL 0

#include <iostream>
#include <vector>
#include <array>
#include <chrono>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

using namespace std;
using namespace boost::random;

const int k = 3;
const int NUM_VERTICES = 10;
int m;

int adj_matrix[NUM_VERTICES][NUM_VERTICES] = {
	{0, 1, 0, 0, 1, 1, 0, 0, 0, 0},
	{1, 0, 1, 0, 0, 0, 1, 0, 0, 0},
	{0, 1, 0, 1, 0, 0, 0, 1, 0, 0},
	{0, 0, 1, 0, 1, 0, 0, 0, 1, 0},
	{1, 0, 0, 1, 0, 0, 0, 0, 0, 1},
	{1, 0, 0, 0, 0, 0, 0, 1, 1, 0},
	{0, 1, 0, 0, 0, 0, 0, 0, 1, 1},
	{0, 0, 1, 0, 0, 1, 0, 0, 0, 1},
	{0, 0, 0, 1, 0, 1, 1, 0, 0, 0},
	{0, 0, 0, 0, 1, 0, 1, 1, 0, 0}
};

const int tenure = 7;
const int num_iterations = 3000;

int s[NUM_VERTICES];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);
uniform_int_distribution<int> random_vertex(0, NUM_VERTICES - 1);

void make_graph(float edge_probability) {  // Populates adj_matrix with a random graph
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
				adj_matrix[i][j] = 1;
				adj_matrix[j][i] = 1;
			}
		}
	}
}


class TabuList {
private:
	int front = 0;
	int ls[tenure][1];
public:
	TabuList() {
		for (int i = 0; i < tenure; i++) {
			ls[i][0] = random_vertex(seed);
			ls[i][1] = random_colour(seed);
		}
	}
	bool contains(int v, int c) {
		for (int i = 0; i < tenure; i++) {
			if (ls[i][0] == v && ls[i][1] == c) {
				return true;
			}
		}
		return false;
	}
	void add(int v, int c) {
		ls[front][0] = v; ls[front][1] = c;
		front = (front + 1) & tenure;
	}
};
TabuList T;

int f(int * x) {
	int num = 0;
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				num++;
			}
		}
	}
	return num;
}

int * A_dict;
int A(int z) {
	if (A_dict[z] != -2) {
		return A_dict[z];
	}
	return z - 1;
}

int neighbour[NUM_VERTICES];
int best[NUM_VERTICES];
void make_move() {  // I've used a rather bold method of resetting neighbour and best (before they're changed) so that we don't need to use copies in the loop
	int best_v = -1;
	int best_c = -1;
	copy(begin(s), end(s), begin(neighbour));
	copy(begin(s), end(s), begin(best));
	for (int v = 0; v < NUM_VERTICES; v++) {
		for (int c = 0; c < k; c++) {
			if(c == 0 && v > 0){
				neighbour[v - 1] = s[v - 1];
			}
			neighbour[v] = c;
			if (f(neighbour) <= A(f(s))) {
				A_dict[f(s)] = f(neighbour) - 1;
				if (f(neighbour) < f(s)) {
					s[v] = c;
					T.add(v, c);
					return;
				}
			} else if (!T.contains(v, c)) {
				if (f(neighbour) < f(s)) {
					s[v] = c;
					T.add(v, c);
					return;
				}
				if (best_v == -1 || f(neighbour) < f(best)) {
					if (best_v != -1) {
						best[best_v] = s[best_v];
					}
					best[v] = c;
					best_v = v;
					best_c = c;
				}
			}
		}
	}
	T.add(best_v, best_c);
	s[best_v] = best_c;
}

void generate_initial_solution() {
	for (int i = 0; i < NUM_VERTICES; i++) {
		s[i] = random_colour(seed);
	}
}

int main(){
	//make_graph(0.5);
	for (int i = 0; i < NUM_VERTICES; i++) {
		for (int j = 0; j < i; j++) {
			m += adj_matrix[i][j];
		}
	}
	generate_initial_solution();
	A_dict = new int[m];
	for (int i = 0; i < m; i++) {
		A_dict[i] = -2;
	}
	int t = 0;
	while (f(s) > 0 && t < num_iterations) {
		make_move();
		t++;
	}
	for (int i = 0; i < NUM_VERTICES; i++) {
		cout << s[i] << " ";
	}
	cout << endl << "Number of conflicts: " << f(s);
	return 0;
}