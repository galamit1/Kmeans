/*** Function Declaration ***/

int get_num_coordinates(char* sinlge_line);
double** init_points (int num_points, int num_coordinates);
double* convert_line_to_point (double* point, char* line);
void free_points_memory (double** points, int num_points);
