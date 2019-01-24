#ifndef KDTREE_CLUSTERING
#define KDTREE_CLUSTERING

//kdtree_clustering
int kdtree(int, int, double*, int, int*, int*, double**, double**);
double* find_centroid(int, double*, int, int);
double* find_variance(int, double*, int, int, double*);
double* calculate_criteria(int, double*, double*);
int bipartition(int, double*, int, int, double*, int[2], int[2]);
double* calculate_boundary(int, double*, int, int);
int assign_values(int*, int[2], int*, int[2], double**, double*[2], double**, double*[2], int*, int*);
double* get_data(double*, int, int);

#endif
