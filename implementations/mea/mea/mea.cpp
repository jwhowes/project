#include <iostream>
#include <algorithm>
#include <boost/random/uniform_real_distribution.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>

using namespace std;
using namespace boost::random;

const int num_vertices = 10;
const int k = 3;

int adj_matrix[num_vertices][num_vertices] = {
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

const int num_iterations = 1;

const int num_agents = 20;
const int q = 500;
const int T_b = 100;
const float r = 1.0f;
const float p = 1.0f;
const int beta = 100;
const int elite_list_size = 5;

int agents[num_agents][num_vertices];
int lifespans[num_agents];
int fitness[num_agents];

int elite_list[elite_list_size][num_vertices];
int elite_fitness[elite_list_size];

mt19937 seed;
uniform_int_distribution<int> random_colour(0, k - 1);
uniform_real_distribution<float> uni(0, 1);
uniform_int_distribution<int> random_elite(0, elite_list_size - 1);

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

void generate_colouring(int * col) {
	for (int v = 0; v < num_vertices; v++) {
		col[v] = random_colour(seed);
	}
}

void generate_agents() {
	for (int i = 0; i < num_agents; i++) {
		generate_colouring(agents[i]);
		lifespans[i] = q;
		fitness[i] = f(agents[i]);
	}
}

void generate_elite_list() {
	for (int i = 0; i < elite_list_size; i++) {
		generate_colouring(elite_list[i]);
		elite_fitness[i] = f(elite_list[i]);
	}
}

void kill_agent(int index) {
	int replace[elite_list_size];
	int replace_num = 0;
	for (int i = 0; i < elite_list_size; i++) {
		if (fitness[index] <= elite_fitness[i]) {
			replace[replace_num] = i;
			replace_num++;
		}
	}
	uniform_int_distribution<int> random_replace(0, replace_num - 1);
	int rep = random_replace(seed);
	copy(begin(agents[index]), end(agents[index]), begin(elite_list[rep]));
	elite_fitness[rep] = fitness[index];
}

void create_agent(int index) {
	int e = random_elite(seed);
	copy(begin(elite_list[e]), end(elite_list[e]), begin(agents[index]));
	lifespans[index] = q;
	fitness[index] = elite_fitness[e];
}

int critical_vertices[num_vertices];
int num_critical;
int tabu_list[num_vertices][num_vertices];
int gamma[num_vertices][num_vertices];
const float lambda = 0.6f;
uniform_int_distribution<int> random_L(0, 9);

void get_critical_vertices(int * col) {
	num_critical = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (gamma[i][col[i]] > 0) {
			critical_vertices[num_critical] = i;
			num_critical++;
		}
	}
}

void update_gamma(int v, int c, int * col) {
	for (int i = 0; i < num_vertices; i++) {
		if (adj_matrix[v][i] == 1) {
			gamma[i][c]++;
			gamma[i][col[v]]--;
		}
	}
}

void tabucol_make_move(int t, int * col) {
	bool initial = true;
	int best_v; int best_c; int best_d;
	get_critical_vertices(col);
	for (int i = 0; i < num_vertices; i++) {
		int v = critical_vertices[i];
		for (int c = 0; c < k; c++) {
			int d = gamma[v][c] - gamma[v][col[v]];
			if (d < 0) {
				update_gamma(v, c, col);
				col[v] = c;
				tabu_list[v][c] = t + random_L(seed) + lambda * num_critical;
				return;
			} else if (tabu_list[v][c] <= t && (initial || d < best_d)) {
				initial = false;
				best_d = d; best_v = v; best_c = c;
			}
		}
	}
	if (!initial) {
		update_gamma(best_v, best_c, col);
		col[best_v] = best_c;
		tabu_list[best_v][best_c] = t + random_L(seed) + lambda * num_critical;
	}
}

int tabucol(int * col, int fitness) {
	int iterations = beta * k * num_vertices / ((k - 1)*fitness);
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < k; j++) {
			tabu_list[i][j] = 0;
			gamma[i][j] = 0;
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < num_vertices; j++) {
			if (adj_matrix[i][j] == 1) {
				gamma[i][col[j]]++;
			}
		}
	}
	for (int t = 0; t < iterations; t++) {
		tabucol_make_move(t, col);
		int t_f = f(col);
		if (t_f == 0) {
			return t_f;
		}
	}
	return f(col);
}

int main(){
	int solution = -1;
	int temp_agent[num_vertices];
	generate_agents();
	for (int t = 0; t < num_iterations; t++) {
		for (int a = 0; a < num_agents; a++) {
			lifespans[a]--;
			if (lifespans[a] <= 0) {
				kill_agent(a);
			}
		}
		if (t % T_b == 0) {
			for (int i = 0; i < num_agents; i++) {
				if (lifespans[i] <= 0) {
					create_agent(i);
					break;
				}
			}
		}
		for (int a = 0; a < num_agents; a++) {
			if (lifespans[a] > 0) {
				copy(begin(agents[a]), end(agents[a]), begin(temp_agent));
				int t_f = tabucol(temp_agent, fitness[a]);
				if (uni(seed) <= r || t_f < fitness[a]) {
					copy(begin(temp_agent), end(temp_agent), begin(agents[a]));
					lifespans[a] += fitness[a] - t_f;
					fitness[a] = t_f;
				}
				if (fitness[a] == 0) {
					solution = a;
				}
			}
		}
		if (solution != -1) {
			break;
		}
	}
	if (solution != -1) {
		for (int v = 0; v < num_vertices; v++) {
			cout << agents[solution][v] << " ";
		}
		cout << endl << "Number of conflicts: " << fitness[solution] << endl;
	} else {
		cout << "No solution found" << endl;
	}
	return 0;
}