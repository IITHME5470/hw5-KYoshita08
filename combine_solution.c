#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Function to read data from a file
void read_solution_file(const char *filename, double *grid_data, int *grid_size) {
    FILE *file = fopen(filename, "r");
    if (!file) {
        fprintf(stderr, "Error: Unable to open file %s\n", filename);
        exit(EXIT_FAILURE);
    }

    int count = 0;
    while (fscanf(file, "%lf", &grid_data[count]) == 1) {
        count++;
    }
    *grid_size = count;
    fclose(file);
}

// Function to write combined data to a single file
void write_combined_solution(const char *output_filename, double *combined_data, int total_size) {
    FILE *file = fopen(output_filename, "w");
    if (!file) {
        fprintf(stderr, "Error: Unable to create output file %s\n", output_filename);
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < total_size; i++) {
        for(int j = 0; j < total_size; j++){
                fprintf(file, "%.6f ", combined_data[i*total_size+j]);
        }
        fprintf(file,"\n");
    }
    fclose(file);
}

int main(int argc, char *argv[]) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <n> <p> <output_file> <solution_file_prefix>\n", argv[0]);
        return EXIT_FAILURE;
    }

    // Read matrix size (n) and number of processors (p)
    int n = atoi(argv[1]);
    int px = atoi(argv[2]);
    int py = atoi(argv[3]);
    const char *output_filename = argv[4];
    const char *file_prefix = argv[5];

    // Calculate total matrix size
    int total_elements = n * n;

    // Allocate memory for combined data
    double *combined_data = malloc(total_elements * sizeof(double));
    if (!combined_data) {
        fprintf(stderr, "Error: Unable to allocate memory for combined data.\n");
        return EXIT_FAILURE;
    }

    int total_size = 0;

    int p = px * py;
    // Loop through each processor's solution file
    for (int i = 0; i < p; i++) {
        char filename[256];
        snprintf(filename, sizeof(filename), "%s_%d_%d_%d.txt", file_prefix, i, px, py); // Generate filename

        double *grid_data = malloc(total_elements / p * sizeof(double));
        if (!grid_data) {
            fprintf(stderr, "Error: Unable to allocate memory for grid data.\n");
            free(combined_data);
            return EXIT_FAILURE;
        }

        int grid_size;
        read_solution_file(filename, grid_data, &grid_size);

        // Append data to combined_data array
        memcpy(&combined_data[total_size], grid_data, grid_size * sizeof(double));
        total_size += grid_size;

        free(grid_data);
    }

    // Write combined data to the output file
    write_combined_solution(output_filename, combined_data, n);

    printf("Combined solution written to %s\n", output_filename);

    free(combined_data);

    return EXIT_SUCCESS;
}