
#define _USE_MATH_DEFINES

#include <iostream>
#include <vector>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <math.h>

using namespace std;
using namespace boost::random;

mt19937 seed;

vector<vector<int>> adj_matrix{
	{0, 1, 1, 0, 1, 0, 0},
	{1, 0, 0, 0, 1, 0, 1},
	{1, 0, 0, 0, 0, 1, 1},
	{0, 0, 0, 0, 0, 0, 1},
	{1, 1, 0, 0, 0, 0, 0},
	{0, 0, 1, 0, 0, 0, 0},
	{0, 1, 1, 1, 0, 0, 0}
};
int n = adj_matrix.size();
int k = 20;

float pa = 0.25;
int num_nests = 50;
int num_iterations = 3000;
bool parasitism_comparison = true;

float beta = 1.5;
float alpha = 1.0;

float sigma_p = pow((tgamma(1+beta)*sin(M_PI*beta/2))/(tgamma((1+beta)/2)*beta*pow(2, (beta-1)/2)), 2 / beta);

vector<vector<int>> nests(num_nests, vector<int>(n, 0));
vector<int> fitness(num_nests, 1);

uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_colour(0, k - 1);
uniform_int_distribution<int> random_vertex(0, n - 1);
normal_distribution<float> normal_q(0, 1);
normal_distribution<float> normal_p(0, sigma_p);

void randomise_colouring(vector<int>& col) {
	for (int i = 0; i < n; i++) {
		col[i] = random_colour(seed);
	}
}

int f(vector<int> col) {
	int num = 0;
	for (int i = 0; i < n; i++) {
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

vector<int> levy_flight(vector<int> col) {
	int M = alpha * levy();
	for (int i = 0; i < M; i++) {
		col[random_vertex(seed)] = random_colour(seed);
	}
	return col;
}

int main(){
	for (int i = 0; i < num_nests; i++) {
		randomise_colouring(nests[i]);
	}
	for (int t = 0; t < num_iterations; t++) {
		for (int i = 0; i < num_nests; i++) {
			vector<int> u_1 = levy_flight(nests[i]);
			int f_u_1 = f(u_1);
			if (f_u_1 <= fitness[i]) {
				nests[i] = u_1;
				fitness[i] = f_u_1;
			}
			if (uni(seed) < pa) {
				vector<int> u_2 = levy_flight(nests[i]);
				int f_u_2 = f(u_2);
				if (!parasitism_comparison || f_u_2 <= fitness[i]) {
					nests[i] = u_2;
					fitness[i] = f_u_2;
				}
			}
		}
	}
	for (int i = 0; i < n; i++) {
		cout << nests[0][i] << " ";
	}
	return 0;
}
