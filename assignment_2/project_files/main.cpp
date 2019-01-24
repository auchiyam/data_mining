#include<iostream>
#include"lsh_cluster.h"
#include"lsh_search.h"

int main() {
	int dim = 2;
	int ndata = 5;
	double* data = new double[dim*ndata]{ 0, 3, 3, 0, -3, 0, 0, -3, 2, 1 };

	int m = 2;
	double** h = new double*[m];
	double w = 3;

	h[0] = new double[2]{-1, 1};
	h[1] = new double[2]{1, 1};

	int* cluster_start;
	int* cluster_size;
	cluster_tree* cluster_hashval;
	int ncluster = -1;

	lsh(dim, ndata, data, m, w, h, &ncluster, &cluster_start, &cluster_size, &cluster_hashval);

	for (int i = 0; i < ncluster; i++) {
		printf("%d\t|\t%d\n", cluster_start[i], cluster_size[i]);
	}

	double* query = new double[2]{0, 4};
	double* result = new double[2];

	lsh_search(dim, ndata, data, m, w, h, ncluster, cluster_start, cluster_size, cluster_hashval, query, result);

	for (int i = 0; i < dim; i++) {
		std::cout << result[i] << ", ";
	}

	std::cout << std::endl;

	return 0;
}
