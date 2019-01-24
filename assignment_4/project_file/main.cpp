#include <iostream>
#include <cmath>
#include <mpi.h>

using namespace std;

int kmean(int, int, double*, int, int*, int*, double*, double**);
int initial_centers(int, int, double*, int, double**);
double* get_point(int, double*, int);
double* get_furthest_point(int, int, double*, double*);
double get_dist(int, double*, double*);
int update_data(int, int, double*, int, int*, int*, double*, double**, double**, int*, int*);
int swap_points(int, double*, int, int);
int search_kmean(int, int, double*, int, int*, int*, double*, double**, double*, double*);
int adjust_data(int, int*, int*, double*, double**);

int main(int argc, char** argv)
{
	MPI_Init(&argc, &argv);

	int dim = 2;
	int ndata = 10;
	double* data = new double[dim*ndata] {1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10};
		
	int k = 3;
		
	int* cluster_start = new int[k];
	int* cluster_size = new int[k];
	double* cluster_radius = new double[k];
	double** centroids = new double*[k];

	int n_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &n_rank);
		
	int count = kmean(dim, ndata, data, k, cluster_start, cluster_size, cluster_radius, centroids);

	if (n_rank == 0) {
		for (int i = 0; i < k; i++) {
			cout << centroids[i][0] << "|" << centroids[i][1] << "|" << cluster_radius[i] << endl;
		}
	
		cout << endl;
		
		for (int i = 0; i < k; i++) {
			cout << cluster_start[i] << "|" << cluster_size[i] << endl;
		}
		
		cout << endl;
	
		for (int i = 0; i < ndata; i++) {
			for (int j = 0; j < dim; j++) {
				cout << data[i*dim+j] << "|";
			}
			cout << endl;
		}
		
		cout << endl;
	}
		
	double* query_pt = new double[dim] {3, 4};
	double* result_pt = new double[dim];
		
	int c = search_kmean(dim, ndata, data, k, cluster_start, cluster_size, cluster_radius, centroids, query_pt, result_pt);

	if (n_rank == 0) {
		cout << c << endl;
		cout << endl;
		
		cout << result_pt[0] << "|" << result_pt[1] << endl;
	}

	MPI_Finalize();

	return 0;
}

int kmean(int dim, int ndata, double *data, int k, int *cluster_start, int *cluster_size, double *cluster_radius, double **centroids) {
	bool no_change = false;
		
	initial_centers(dim, ndata, data, k, centroids);

	int n_rank;

	MPI_Comm_rank(MPI_COMM_WORLD, &n_rank);

	if (n_rank == 0) {
		int* cluster_assign = new int[ndata];
		for (int i = 0; i < ndata; i++) {
			cluster_assign[i] = -1;
		}
		
		while (!no_change) {
			//initialization
			no_change = true;
					
			double **sum = new double*[k];
			for (int i = 0; i < k; i++) {
				sum[i] = new double[k];
				for (int j = 0; j < dim; j++) {
					sum[i][j] = 0;
				}
			}
			
			int *count = new int[k];
			for (int i = 0; i < k; i++) {
				count[i] = 0;
			}
				
			//assign all data points to a cluster
			for (int i = 0; i < ndata; i++) {
				double* point = get_point(dim, data, i);
			
				double min_dist = -1;
				int cluster = -1;
				//find cluster that's closest to point i
				for (int j = 0; j < k; j++) {
					double dist = get_dist(dim, point, centroids[j]);
						
					if (min_dist == -1 || dist < min_dist) {
						cluster = j;
						min_dist = dist;
					}
				}
					
				//check if the closest cluster is the same
				if (cluster != cluster_assign[i]) {
					no_change = false;
						
					cluster_assign[i] = cluster;
				}
					
				//add the point to the sum and increase counter to find the average later
				for (int j = 0; j < dim; j++) {
					sum[cluster][j] += point[j];
				}
				count[cluster]++;
			}
		
			//update everything after every iteration
			update_data(dim, ndata, data, k, cluster_start, cluster_size, cluster_radius, centroids, sum, count, cluster_assign);
		}

		int count = adjust_data(k, cluster_start, cluster_size, cluster_radius, centroids);

		MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(cluster_start, k, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(cluster_size, k, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(cluster_radius, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int i = 0; i < k; i++) {
			MPI_Bcast(centroids[i], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		return count;
	} else {
		int count;

		MPI_Bcast(&count, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(cluster_start, k, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(cluster_size, k, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(cluster_radius, k, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		for (int i = 0; i < k; i++) {
			centroids[i] = new double[dim];
			MPI_Bcast(centroids[i], dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		}

		return count;
	}
}

int initial_centers(int dim, int ndata, double *data, int k, double **cluster_centroid) {
	cluster_centroid[0] = get_point(dim, data, rand() % ndata);
	cluster_centroid[1] = get_furthest_point(dim, ndata, data, cluster_centroid[0]);
		
	for (int i = 2; i < k; i++) {
		double max_dist = -1;
		double* max_point;
		double* reduced_max_point;

		int n_tasks, n_rank;

		MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
		MPI_Comm_rank(MPI_COMM_WORLD, &n_rank);

		int iterate_count = ndata / n_tasks; //assumption: ndata can be evenly divided by amount of parallel tasks
				
		//find maximum minimum distance in the subset of nodes
		for (int j = iterate_count * n_rank; j < iterate_count * (n_rank+1); j++) {
			//distance for closest
			double min_dist = -1;
			double* min_point;
						
			for (int l = 0; l < i; l++) {
				double* point = get_point(dim, data, j);
				double dist = get_dist(dim, cluster_centroid[l], point);

				//update shortest distance for point j
				if (min_dist == -1 || dist < min_dist) {
					min_dist = dist;
					min_point = point;
				}
			}
						
			if (min_dist > max_dist) {
				max_dist = min_dist;
				max_point = min_point;
			}
		}

		if (n_rank == 0) {
			double* dist = new double[n_tasks];
			double** point = new double*[n_tasks];

			//wait for all the ranks to be done
			for (int j = 1; j < n_tasks; j++) {
				MPI_Status status;
				MPI_Recv(&dist[j], 1, MPI_DOUBLE, j, 0, MPI_COMM_WORLD, &status);

				point[j] = new double[dim];
				MPI_Recv(point[j], dim, MPI_DOUBLE, j, 1, MPI_COMM_WORLD, &status);
			}

			//find the max dist out of them all
			for (int j = 1; j < n_tasks; j++) {
				if (dist[j] > max_dist) {
					max_dist = dist[j];
					max_point = point[j];
				}
			}

			cluster_centroid[i] = max_point;
		} else {
			MPI_Send(&max_dist, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(max_point, dim, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
		}
	}

	return 0;
}

double* get_point(int dim, double* data, int loc) {
	double* pt = new double[dim];
	for (int i = 0; i < dim; i++) {
		pt[i] = data[loc*dim+i];
	}
		
	return pt;
}

double* get_furthest_point(int dim, int ndata, double* data, double* point) {
	double max_dist = -1;
	double* max_point;
	for (int i = 0; i < ndata; i++) {
		double* pt = get_point(dim, data, i);
		
		double dist = get_dist(dim, point, pt);
		
		if (dist > max_dist) {
			max_dist = dist;
			max_point = pt;
		}
	}
		
	return max_point;
}

double get_dist(int dim, double* pt_1, double* pt_2) {
	double sum = 0;
		
	for (int i = 0; i < dim; i++) {
		sum += pow((pt_1[i] - pt_2[i]), 2);
	}
		
	return sqrt(sum);
}

int update_data(int dim, int ndata, double* data, int k, int* cluster_start, int* cluster_size, double* cluster_radius, double** centroid, double** sum_points, int* count, int* cluster_assign) {
	//initialization
	int* counter = new int[k];
	int* ind = new int[k];
	for (int i = 0; i < k; i++) {
		counter[i] = 0;
		ind[i] = 0;
	}
		
	//count the size of each clusters
	for (int i = 0; i < ndata; i++) {
		counter[cluster_assign[i]]++;
	}
		
	int sum = 0;
	//update cluster start and size
	for (int i = 0; i < k; i++) {
		cluster_start[i] = sum;
		cluster_size[i] = counter[i];
			
		sum += counter[i];
	}
	
	int curr_cluster = 0;
	
	//update data
	for (int i = 0; i < ndata; i++) {
		bool proceed = false;
		//check if the current i is on the next cluster or not
		if (i == cluster_start[curr_cluster+1]) {
			curr_cluster++;
				
			//skip all the slots we know is already at the proper area
			i += ind[curr_cluster];
			
			//make sure this jump didn't go to the next cluster
			if (i >= cluster_start[curr_cluster]) {
				i--;
				continue;
			}
		}
		
		//keep swapping at i until point at i is correct cluster
		while (!proceed) {
			int cluster = cluster_assign[i];
	
			//point is in wrong place
			if (cluster != curr_cluster) {
				int correct_index = cluster_start[cluster] + ind[cluster];
				swap_points(dim, data, i, correct_index);
				
				cluster_assign[i] = cluster_assign[correct_index];
				cluster_assign[correct_index] = cluster;
				
				ind[cluster]++;
			}
			//the point is in the right place
			else {
				proceed = true;
			}
		}
	}
	
	//update centroids
	for (int i = 0; i < k; i++) {
		for (int j = 0; j < dim; j++) {
			centroid[i][j] = sum_points[i][j] / count[i];
		}
	}
	
	//update radius
	for (int i = 0; i < k; i++) {
		double max_dist = -1;
		double* max_point;
		for (int j = cluster_start[i]; j < cluster_size[i] + cluster_start[i]; j++) {
			double* centroid_pt = centroid[i];
			double* pt = get_point(dim, data, j);
			
			double dist = get_dist(dim, centroid_pt, pt);
			
			if (dist > max_dist) {
				max_dist = dist;
				max_point = pt;
			}
		}
			
		cluster_radius[i] = max_dist;
	}
	
	return 0;
}

int adjust_data(int k, int* cluster_start, int* cluster_size, double* cluster_radius, double** centroids) {
	int empty_count = 0;
	for (int i = 0; i+empty_count < k; i++) {
		if (cluster_size[i] == 0) {
			empty_count++;
			
			if (i+empty_count == k) {
				break;
			}
		}
		
		//remove all empty clusters and shift the rest to fill the gap
		cluster_start[i] = cluster_start[i+empty_count];
		cluster_size[i] = cluster_size[i+empty_count];
		cluster_radius[i] = cluster_radius[i+empty_count];
		centroids[i] = centroids[i+empty_count];
	}
		
	return k - empty_count; // new cluster count
}

int swap_points(int dim, double* data, int i_1, int i_2) {
	for (int i = 0; i < dim; i++) {
		double tmp = data[i_1*dim+i];
		data[i_1*dim+i] = data[i_2*dim+i];
		data[i_2*dim+i] = tmp;
	}
	return 0;
}

int search_kmean(int dim, int ndata, double *data, int k, int *cluster_start, int *cluster_size, double *cluster_radius, double **centroids, double *query_pt, double *result_pt) {
	int* viable_clusters = new int[k];
	int min_cluster = -1;
	double min_dist = -1;
	int cc = 0;
	int n_tasks, n_rank;
	int total_checked = 0;
	int checked = 0;

	MPI_Comm_size(MPI_COMM_WORLD, &n_tasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &n_rank);

	//find all the clusters query point can be part of
	if (n_rank < k) {
		for (int i = n_rank; i < k; i += n_tasks) {
			double dist = get_dist(dim, centroids[i], query_pt);
			
			//if the dist between centroid and query is smaller than the maximum radius, then the cluster is relevant
			if (dist < cluster_radius[i]) {
				viable_clusters[cc] = i;
				cc++;
			}
			
			if (min_dist == -1 || dist < min_dist) {
				min_dist = dist;
				min_cluster = i;
			}
		}

		//the point is not in any of the clusters, so choose the closest cluster instead
		if (cc == 0) {
			cc++;
	
			viable_clusters[0] = min_cluster;
		}

		min_dist = -1;
		double* min_pt;
	
		//look at every points in the relevant clusters and get the closest points
		for (int i = 0; i < cc; i++) {
			int cluster = viable_clusters[i];
		
			for (int j = cluster_start[cluster]; j < cluster_start[cluster] + cluster_size[cluster]; j++) {
				double* pt = get_point(dim, data, j);
				double dist = get_dist(dim, query_pt, pt);
				
				if (min_dist == -1 || dist < min_dist) {
					min_dist = dist;
					min_pt = pt;
				}
				
				checked++;
			}
		}

		if (n_rank == 0) {
			int count = n_tasks;
			if (k <= n_tasks) {
				count = k;
			}
				
			for (int i = 1; i < count; i++) {
				double* pt = new double[dim];
				double dist;

				MPI_Recv(pt, dim, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
				MPI_Recv(&dist, 1, MPI_DOUBLE, i, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

				if (dist < min_dist) {
					min_dist = dist;
					min_pt = pt;
				}
			}

			//set result_pt to the point found
			for (int i = 0; i < dim; i++) {
				result_pt[i] = min_pt[i];
			}

			MPI_Bcast(result_pt, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Allreduce(&checked, &total_checked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		} else {
			MPI_Send(min_pt, dim, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
			MPI_Send(&min_dist, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);

			MPI_Bcast(result_pt, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
			MPI_Allreduce(&checked, &total_checked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		}
	} else {
		MPI_Bcast(result_pt, dim, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		MPI_Allreduce(&checked, &total_checked, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}

	return total_checked;
}
