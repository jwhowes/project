#include <iostream>

const int num_vertices = 10;
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
int adj_list[num_vertices][num_vertices] = {
	{1, 4, 5, 0, 0, 0, 0, 0, 0, 0},
	{0, 2, 6, 0, 0, 0, 0, 0, 0, 0},
	{1, 3, 7, 0, 0, 0, 0, 0, 0, 0},
	{2, 4, 8, 0, 0, 0, 0, 0, 0, 0},
	{0, 3, 9, 0, 0, 0, 0, 0, 0, 0},
	{0, 7, 8, 0, 0, 0, 0, 0, 0, 0},
	{1, 8, 9, 0, 0, 0, 0, 0, 0, 0},
	{2, 5, 9, 0, 0, 0, 0, 0, 0, 0},
	{3, 5, 6, 0, 0, 0, 0, 0, 0, 0},
	{4, 6, 7, 0, 0, 0, 0, 0, 0, 0}
};
int adj_list_length[num_vertices] = { 3, 3, 3, 3, 3, 3, 3, 3, 3, 3 };

int tabu_list[num_vertices][2];

const int num_iterations = 10;

void inline generate_initial_solution(int orientation[num_vertices][num_vertices]) {
	for (int i = 0; i < num_vertices; i++) {
		for(int j = i; j < num_vertices; j++){
			if (adj_matrix[i][j] == 1) {
				// j > i (hence orientation i <- j)
				orientation[i][j] = -1;
				orientation[j][i] = 1;
			} else {
				orientation[i][j] = orientation[i][j] = 0;
			}
		}
	}
}

int populate_distances(int orientation[num_vertices][num_vertices], int * d_plus, int * d_minus) {
	int indegree[num_vertices];
	int outdegree[num_vertices];
	int queue[num_vertices];
	int queue_length;
	int lambda = 0;
	for (int i = 0; i < num_vertices; i++) {
		d_minus[i] = 0;
		d_plus[i] = 0;
		indegree[i] = 0;
		outdegree[i] = 0;
	}
	queue_length = 0;
	for (int i = 0; i < num_vertices; i++) {
		for (int j = 0; j < adj_list_length[i]; j++) {
			if (orientation[i][adj_list[i][j]] == -1) {
				indegree[i]++;
			} else {
				outdegree[i]++;
			}
		}
	}
	for (int i = 0; i < num_vertices; i++) {
		if (indegree[i] == 0) {
			queue[queue_length] = i;
			queue_length++;
		}
	}
	int num_checked = 0;
	while (num_checked < num_vertices) {
		int i = queue[num_checked];
		for (int j = 0; j < adj_list_length[i]; j++) {
			if (orientation[i][adj_list[i][j]] == 1) {
				indegree[j] -= 1;
				if (d_minus[i] + 1 > d_minus[j]) {
					d_minus[j] = d_minus[i] + 1;
					if (d_minus[i] > lambda) {
						lambda = d_minus[i];
					}
				}
				if (indegree[j] == 0) {
					queue[queue_length] = j;
					queue_length++;
				}
			}
		}
		num_checked++;
	}
	queue_length = 0;
	for (int i = 0; i < num_vertices; i++) {
		if (outdegree[i] == 0) {
			queue[queue_length] = i;
			queue_length++;
		}
	}
	num_checked = 0;
	while (num_checked < num_vertices) {
		int i = queue[num_checked];
		for (int j = 0; j < adj_list_length[i]; j++) {
			if (orientation[i][adj_list[i][j]] == -1) {
				outdegree[j] -= 1;
				if (d_plus[i] + 1 > d_plus[j]) {
					d_plus[j] = d_plus[i] + 1;
				}
				if (outdegree[j] == 0) {
					queue[queue_length] = j;
					queue_length++;
				}
			}
		}
		num_checked++;
	}
	return lambda;
}

int main(){
	bool found[num_vertices];
	int on_longest_path[num_vertices];
	int on_longest_length;
	int orientation[num_vertices][num_vertices];
	int d_plus[num_vertices];
	int d_minus[num_vertices];
	generate_initial_solution(orientation);
	int indegree[num_vertices];
	int queue[num_vertices];
	int queue_length;
	for (int i = 0; i < num_vertices; i++) {
		tabu_list[i][0] = 0;
		tabu_list[i][1] = 0;
	}
	int lambda = populate_distances(orientation, d_plus, d_minus);
	int best_lambda = lambda;
	for (int t = 0; t < num_iterations; t++) {
		on_longest_length = 0;
		for (int i = 0; i < num_vertices; i++) {
			found[i] = false;
		}
		for (int i = 0; i < num_vertices; i++) {
			for (int j = i; j < num_vertices; j++) {
				if ((orientation[i][j] == 1 && d_minus[i] + d_plus[j] == lambda) || (orientation[i][j] == -1 && d_minus[j] + d_plus[i] == lambda)) {
					if (!found[i]) {
						found[i] = true;
						on_longest_path[on_longest_length] = i;
						on_longest_length++;
					}
					if (!found[j]) {
						found[j] = true;
						on_longest_path[on_longest_length] = j;
						on_longest_length++;
					}
				}
			}
		}
	}
	return 0;
}