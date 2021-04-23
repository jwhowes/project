#define _SECURE_SCL 0

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <array>
#include <math.h>
#include <chrono>
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace boost::random;

const string graph_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/graphs/";
const string results_directory = "C:/Users/taydo/OneDrive/Documents/computer_science/year3/project/implementations/results/";

const int num_vertices = 500;
int adj_matrix[num_vertices][num_vertices];
int adj_list[num_vertices][num_vertices];
int adj_list_length[num_vertices];

int k;

const int num_iterations = 3000;
const auto duration = chrono::minutes{ 5 };


const int num_cols = 50;

int cols[num_cols][num_vertices];
int fitness[num_cols];

struct Reference {
	int col[num_vertices];
	int fitness;
	int diversity;
};

const int b = 5;
const int d = 5;
Reference ref_set[b + d];
int worst_best_set = 0;
int most_diverse_div_set = 0;

int best_colouring[num_vertices];

mt19937 seed;
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
uniform_int_distribution<int> random_colour;
uniform_int_distribution<int> random_ref_set(0, b + d - 1);

void get_random_colouring(int * col) {
	for (int i = 0; i < num_vertices; i++) {
		col[i] = random_colour(seed);
	}
}

bool valid(int v, int c, int * col) {  // Returns whether or not vertex v can be coloured colour c in colouring col (legally)
	for (int i = 0; i < adj_list_length[v]; i++) {
		if (col[adj_list[v][i]] == c) {  // If v is adjacent to some vertex i coloured c then this is not a valid assignment
			return false;
		}
	}
	return true;  // If no such i can be found then the assignment is valid
}

int chromatic_bound() {
	int ret = 0;
	for (int i = 0; i < num_vertices; i++) {
		int c = 0;
		while (true) {
			if (valid(i, c, best_colouring)) {
				best_colouring[i] = c;
				if (c > ret) {
					ret = c;
				}
				break;
			}
			c++;
		}
	}
	return ret + 1;
}

int f(int * x) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && x[i] == x[j]) {
				num++;
			}
		}
	}
	return num;
}

int main(){
	cout << "SS\n";
	k = chromatic_bound();
	random_colour = uniform_int_distribution<int>(0, k - 1);
	for (int i = 0; i < num_cols; i++) {
		get_random_colouring(cols[i]);
		fitness[i] = f(cols[i]);
	}
	return 0;
}