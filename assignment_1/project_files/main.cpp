#include <iostream>
#include "kdtree_clustering.h"
#include "kdtree_search.h"

int main() {
	double* data = new double[20] {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 11};
	int ndata = 10;
	int dim = 2;

	for (int i = 0; i < ndata * dim; i++) {
		data[i] = rand() % 100;
	}

	int kk = 2;
	int* cluster_start = new int[kk];
	int* cluster_size = new int[kk];
	double** cluster_bdry = new double*[kk];
	double** centroid = new double*[kk];

	kdtree(dim, ndata, data, 2, cluster_start, cluster_size, cluster_bdry, centroid);

	std::cout << "Every points:" << std::endl;
	for (int i = 0; i < ndata; i++) {
		for (int j = 0; j < dim; j++) {
			std::cout << data[i*dim+j] << "|";
		}	
		std::cout << std::endl;
	}

	for (int i = 0; i < kk; i++) {
		printf("%d: %d | size: %d\n", i, cluster_start[i], cluster_size[i]);
	}

	std::cout << std::endl;

	for (int i = 0; i < kk; i++) {
		for (int j = 0; j < 2 * dim; j += 2) {
			printf("%d: %f < dim(%d) < %f\n", i, cluster_bdry[i][j], j/2, cluster_bdry[i][j + 1]);
		}
	}

	std::cout << std::endl;

	double* query_pt = new double[dim] { 7, 3 };
	double* result_pt = new double[dim];

	int count = kdsearch(dim, ndata, data, kk, cluster_start, cluster_size, cluster_bdry, query_pt, result_pt);

	printf("query:\t(%f,%f)\n", query_pt[0], query_pt[1]);
	printf("result:\t(%f,%f)\n", result_pt[0], result_pt[1]);
	printf("count:\t%d", count);

	std::cout << std::endl;

	return 0;
}
