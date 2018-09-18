#include<iostream>

int kdtree(int, int, double*, int, int*, int*, double**, double**);
int* find_centroid(int, int, double*, int, int);
int* find_variance(int, int, double*, int, int, int*);
int* calculate_boundary(int, int, double*, int, int, int*, int*);
int bipartition(int, int, int, double*, int[2], int[2], double*[2], double*[2]);

int kdtree(int dim, int ndata, double* data, int kk, int* cluster_start, int* cluster_size, double** cluster_bdry, double** centroid) {
	//initialize all the outputs
	cluster_start = new int[kk];
	cluster_size = new int[kk];
	cluster_bdry = new int*[kk];
	centroid = new double*[kk];

	for (int i = 0; i < kk; i++) {
		cluster_bdry[i] = new int[2 * dim];
		centroid[i] = new double[dim]
	}

	//initialize all the other variables
	bool done = false;

	//loop while we're not done
	while(!done) {
		int* curr_centroid = find_centroid(dim, ndata, data, cluster_start, cluster_size);
		int* curr_variance = find_variance(dim, ndata, data, cluster_start, cluster_size, curr_centroid);
		int* curr_boundary = calculate_boundary(dim, ndata, data, cluster_start, cluster_size, curr_variance, curr_centroid);
		int* new_cluster_start = new int[2];
		int* new_cluster_size = new int[2];
		int* new_cluster_bdry = new int*[2];
		bipartition(dim, , , data, new_cluster_start, new_cluster_size, new_cluster_bdry
	}
}


int* find_centroid(int dim, int ndata, double* data, int cluster_start, int cluster_size) {

}

int* find_variance(int dim, int ndata, double*data, int cluster_start, int cluster_size) {

}

int* calculate_boundary(int dim, int ndata, double* data, int cluster_start, int cluster_size, int* variance,  int* centroid) {

}

int bipartition(int dim, int i0, int in, double* data, int cluster_start[2], int cluster_size[2], double* cluster_bdry[2], double* centroid[2]) {

}
