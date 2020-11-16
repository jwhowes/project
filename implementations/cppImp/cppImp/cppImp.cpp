
#define _USE_MATH_DEFINES
#define _SECURE_SCL 0
#define num_vertices 100
#define num_nests 50

#include <iostream>
//#include <array>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <chrono>
#include <math.h>

using namespace std;
using namespace boost::random;

mt19937 seed;

/*int adj_matrix[num_vertices][num_vertices] = {
	{0, 1, 1, 0, 1, 0, 0},
	{1, 0, 0, 0, 1, 0, 1},
	{1, 0, 0, 0, 0, 1, 1},
	{0, 0, 0, 0, 0, 0, 1},
	{1, 1, 0, 0, 0, 0, 0},
	{0, 0, 1, 0, 0, 0, 0},
	{0, 1, 1, 1, 0, 0, 0}
};*/
int adj_matrix[num_vertices][num_vertices];
const int k = 3;

const float pa = 0.25;
const int num_iterations = 3000;
const bool parasitism_comparison = true;

const float beta = 1.5;
const float alpha = 1.0;

const float sigma_p = pow((tgamma(1+beta)*sin(M_PI*beta/2))/(tgamma((1+beta)/2)*beta*pow(2, (beta-1)/2)), 2 / beta);

int nests[num_nests][num_vertices];
int fitness[num_nests];

uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);
uniform_int_distribution<int> random_vertex(0, num_vertices - 1);
normal_distribution<float> normal_q(0, 1);
normal_distribution<float> normal_p(0, sigma_p);

void get_graph(float edge_probability) {
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (uni(seed) < edge_probability) {
				adj_matrix[i][j] = 1;
				adj_matrix[j][i] = 1;
			}
		}
	}
}

void randomise_colouring(int * col) {
	for (int i = 0; i < num_vertices; i++) {
		col[i] = random_colour(seed);
	}
}

int f(int * col) {
	int num = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < i; j++) {
			if (adj_matrix[i][j] == 1 && col[i] == col[j]) {
				num++;
			}
		}
	}
	return num;
}

float levy() {
	float p = normal_p(seed);
	float q = normal_q(seed);
	return p / pow(abs(q), 1 / beta);
}

void levy_flight(int * col) {
	int M = alpha * levy();
	for (int i = 0; i < M; i++) {
		col[random_vertex(seed)] = random_colour(seed);
	}
}

void set(int * nest, int * replacement) {
	for(int i = 0; i < num_vertices; i++){
		nest[i] = replacement[i];
	}
}

int main(){
	get_graph(0.5);
	int u_1[num_vertices];
	int u_2[num_vertices];
	auto start = chrono::high_resolution_clock::now();
	for (int i = 0; i < num_nests; i++) {
		randomise_colouring(nests[i]);
	}
	for (int t = 0; t < num_iterations; t++) {
		for (int i = 0; i < num_nests; i++) {
			copy(begin(nests[i]), end(nests[i]), begin(u_1));
			levy_flight(u_1);
			int f_u_1 = f(u_1);
			if (f_u_1 <= fitness[i]) {
				set(nests[i], u_1);
				fitness[i] = f_u_1;
			}
			if (uni(seed) < pa) {
				copy(begin(nests[i]), end(nests[i]), begin(u_2));
				levy_flight(u_2);
				int f_u_2 = f(u_2);
				if (!parasitism_comparison || f_u_2 <= fitness[i]) {
					set(nests[i], u_2);
					fitness[i] = f_u_2;
				}
			}
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		cout << nests[0][i] << " ";
	}
	cout << endl << "Time taken: " << chrono::duration_cast<chrono::microseconds>(chrono::high_resolution_clock::now() - start).count();
	return 0;
}
