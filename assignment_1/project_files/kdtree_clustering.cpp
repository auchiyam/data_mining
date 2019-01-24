#include <iostream>
#include <cmath>
#include "kdtree_clustering.h"
int kdtree(int dim, int ndata, double* data, int kk, int* cluster_start, int* cluster_size, double** cluster_bdry, double** centroid) {
	//initialize all the other variables
	bool done = false;
	int empty_slot = 1;
	int curr_cluster = 0;
	int total_cluster = 1;

	centroid[0] = find_centroid(dim, data, 0, ndata);
	cluster_start[0] = 0;
	cluster_size[0] = ndata;

	//loop while we're not done
	while(!done) {
		double* curr_centroid = centroid[curr_cluster];
		double* curr_variance = find_variance(dim, data, cluster_start[curr_cluster], cluster_size[curr_cluster], curr_centroid);
		double* curr_criteria = calculate_criteria(dim, curr_variance, curr_centroid);

		int* new_cluster_start = new int[2];
		int* new_cluster_size = new int[2];
		double** new_cluster_bdry = new double*[2];
		double** new_centroids = new double*[2];

		bipartition(dim, data, cluster_start[curr_cluster], cluster_size[curr_cluster], curr_criteria, new_cluster_start, new_cluster_size);
		new_cluster_bdry[0] = calculate_boundary(dim, data, new_cluster_start[0], new_cluster_size[0]);
		new_cluster_bdry[1] = calculate_boundary(dim, data, new_cluster_start[1], new_cluster_size[1]);
		new_centroids[0] = find_centroid(dim, data, new_cluster_start[0], new_cluster_size[0]);
		new_centroids[1] = find_centroid(dim, data, new_cluster_start[1], new_cluster_size[1]);

		assign_values(cluster_start, new_cluster_start, cluster_size, new_cluster_size, cluster_bdry, new_cluster_bdry, centroid, new_centroids, &curr_cluster, &empty_slot);

		if (curr_cluster == total_cluster) {
			total_cluster *= 2;
		       	curr_cluster = 0;
		}

		if (total_cluster >= kk) {
			done = true;
		}
	}

	return 0;
}


double* find_centroid(int dim, double* data, int cluster_start, int cluster_size) {
	double* average = new double[dim];

	for (int i = cluster_start; i < cluster_start + cluster_size; i++) {
		for (int j = 0; j < dim; j++) {
			//make sure average is initialized on the first run
			if (i == 0) {
				average[j] = 0;
			}
			average[j] += data[i * dim + j];
		}
	}

	for (int i = 0; i < dim; i++) {
		average[i] /= cluster_size;
	}

	return average;
}

double* find_variance(int dim, double* data, int cluster_start, int cluster_size, double* centroid) {
	double* average = new double[dim];

	for (int i = cluster_start; i < cluster_start + cluster_size; i++) {
		for (int j = 0; j < dim; j++) {
			//make sure average is initialized on the first run
			if (i == 0) {
				average[j] = 0;
			}
			average[j] += std::pow((data[i * dim + j] - centroid[j]), 2);
		}
	}

	for (int i = 0; i < dim; i++) {
		average[i] /= cluster_size;
	}

	return average;
}

double* calculate_criteria(int dim, double* variance,  double* centroid) {
	int max_dim = -1;
	double max = 0;

	//get the dimension with highest number
	for (int i = 0; i < dim; i++) {
		if (variance[i] > max) {
			max = variance[i];
			max_dim = i;
		}
	}

	return new double[2] { (double)max_dim, centroid[max_dim] };
}

int bipartition(int dim, double* data, int i0, int in, double* criteria, int cluster_start[2], int cluster_size[2]) {
	double* cluster_2 = new double[in * dim];
	int iterator_1 = 0;
	int iterator_2 = 0;

	//Separate the points into two groups
	for (int i = i0; i < i0 + in; i++) {
		double* curr_point = get_data(data, dim, i);

		if (curr_point[(int)criteria[0]] < criteria[1]) {
			for (int j = 0; j < dim; j++) {
				data[(i0 + iterator_1) * dim + j] = curr_point[j];
			}

			iterator_1++;
		}
		else {
			for (int j = 0; j < dim; j++) {
				cluster_2[iterator_2 * dim + j] = curr_point[j];
			}

			iterator_2++;
		}
	}

	//Add the second cluster to the end of data
	for (int i = 0; i < iterator_2 * dim; i++) {
		data[(i0 + iterator_1) * dim + i] = cluster_2[i];
	}

	cluster_start[0] = i0;
	cluster_start[1] = i0 + iterator_1;

	cluster_size[0] = iterator_1;
	cluster_size[1] = iterator_2;
	//success

	return 0;
}

double* calculate_boundary(int dim, double* data, int cluster_start, int cluster_size)
{
	double* min = new double[dim];
	double* max = new double[dim];

	double* new_boundary = new double[2 * dim];

	for (int i = 0; i < dim; i++) {
		min[i] = data[i + cluster_start * dim];
		max[i] = data[i + cluster_start * dim];
	}

	for (int i = cluster_start; i < cluster_start + cluster_size; i++) {
		for (int j = 0; j < dim; j++) {
			double val = data[i * dim + j];
			if (val < min[j]) {
				min[j] = val;
			}

			if (val > max[j]) {
				max[j] = val;
			}
		}
	}

	for (int i = 0; i < 2 * dim; i += 2) {
		new_boundary[i] = min[i/2];
		new_boundary[i+1] = max[i/2];
	}	

	return new_boundary;
}

double* get_data(double* data, int dim, int i) {
	double* data_point = new double[dim];

	for (int j = 0; j < dim; j++)
	{
		data_point[j] = data[i*dim+j];
	}

	return data_point;
}

int assign_values(int* cluster_start, int new_cluster_start[2], int* cluster_size, int new_cluster_size[2], double** cluster_bdry, double* new_cluster_bdry[2], double** centroid, double* new_centroids[2], int* curr_cluster, int* empty_slot) {
	cluster_start[*curr_cluster] = new_cluster_start[0];
	cluster_start[*empty_slot] = new_cluster_start[1];

	cluster_size[*curr_cluster] = new_cluster_size[0];
	cluster_size[*empty_slot] = new_cluster_size[1];

	cluster_bdry[*curr_cluster] = new_cluster_bdry[0];
	cluster_bdry[*empty_slot] = new_cluster_bdry[1];

	centroid[*curr_cluster] = new_centroids[0];
	centroid[*empty_slot] = new_centroids[1];

	*curr_cluster += 1;
	*empty_slot += 1;

	return 0;
}
