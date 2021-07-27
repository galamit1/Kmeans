/*** Function Declaration ***/

typedef struct Cluster {
    int cluster_size;
    int cluster_index;
    double* prev_centroids;
    double* curr_centroids;
    double* sum;
} Cluster;

int get_num_coordinates(char* sinlge_line);
double* convert_line_to_point (double* point, char* line);
Cluster* make_cluster (const double* point, const int num_coordinates, const int index);
Cluster** init_k_clusters (double** points, int k, int num_coordinates);
void add_point_to_cluster (Cluster* cluster, const double* point, int num_coordinates);
void find_minimum_distance (Cluster** clusters, double* point, int k, int num_coordinates);
double update_cluster_centroids_and_sum (Cluster* cluster, int num_coordinates);
int update_all_clusters (Cluster** clusters, int k, int num_coordinates);
void free_clusters_memory (Cluster** clusters, int k);
double get_distance (Cluster* cluster, const double* point, int num_coordinates);
