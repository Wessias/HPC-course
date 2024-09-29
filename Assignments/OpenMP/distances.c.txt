#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <string.h>

#define MAX_BLOCK_SIZE 200000  // Number of cells in one block. Currently using 2 blocks so 
			       // block size has to be less than 5242800/24 = 218450 
#define MAX_DISTANCES  5000   // Max distance is approx 37.00 so possible values after round
			      // is for sure less than 5000.

typedef struct {
	float x, y, z;  
} Cell;

int distance_count[MAX_DISTANCES] = {0};
// Used to store frequency of distances. (GLOBAL VARIABLE)
// Indices represent the lengths i.e. index 1234
// stores the count of how often 12.34 has shown up.

/*
 * Calculates distance between two 3D points in the usual manner.
 * Truncates the resulting float with 2 decimals by multiplying by 100 and casting to short.
 * Uses the truncated length to increase the local_distance_count.
 */
static inline void calculate_distance(Cell c1, Cell c2, int *local_distance_count) {
	float dx = c1.x - c2.x;
	float dy = c1.y - c2.y;
	float dz = c1.z - c2.z;
	short rounded_dist = (int)(sqrtf(dx * dx + dy * dy + dz * dz) * 100);  // Use sqrtf for float
	local_distance_count[rounded_dist]++;

}

/*
 * Calculates distances between all unique combinations of cells in a block.
 * Can skip a lot of combinations that will occur in the processing of two different blocks.
 * Parallel outer loop and vectorized inner loop.
 * Automatically updates distance_count.
 */
static inline void process_block(Cell *block, int num_cells) {
#pragma omp parallel
	{
		int local_distance_count[MAX_DISTANCES] = {0};
#pragma omp for schedule(dynamic, 100)  
		for (int i = 0; i < num_cells; ++i) {
#pragma omp simd
			for (int j = i + 1; j < num_cells; ++j) {
				calculate_distance(block[i], block[j], local_distance_count);
			}
		}

#pragma omp critical
		{
			for (int k = 0; k < MAX_DISTANCES; ++k){
				distance_count[k] += local_distance_count[k];
			}

		}
	}
}

/*
 * Calculates distances between two different blocks.
 * Usual looping of two "matrices".
 * Parallel outer loop and vectorized inner.
 * Automatically updates distance_count.
 */
static inline void process_inter_block(Cell *main_block, int num_current, Cell *secondary_block, int num_secondary) {
#pragma omp parallel
	{


		int local_distance_count[MAX_DISTANCES] = {0};
#pragma omp for schedule(dynamic, 100)
		for (int i = 0; i < num_current; ++i) {
#pragma omp simd
			for (int j = 0; j < num_secondary; ++j) {
				calculate_distance(main_block[i], secondary_block[j], local_distance_count);

			}
		}
#pragma omp critical
		{
			for (int k = 0; k < MAX_DISTANCES; ++k){
				distance_count[k] += local_distance_count[k];
			}

		}
	}
}

static inline void output_distances(int *distance_count) {
	for (int i = 1; i < MAX_DISTANCES; ++i) {
		if (distance_count[i] > 0) {
			printf("%05.2f %d\n", i / 100.0f, distance_count[i]);  // Use 100.0f for float division
		}
	}
}

int main(int argc, char *argv[]) {
	int num_threads = 1;
	if (argc != 2) {
		num_threads = 1;
	} else {
		num_threads = atoi(argv[1] + 2);  	}

	omp_set_num_threads(num_threads);     

	FILE *file = fopen("cells", "r");
	if (!file) {
		perror("Error opening file");
		return 1;
	}


	Cell *main_block = (Cell *)malloc(MAX_BLOCK_SIZE * sizeof(Cell));
	Cell *secondary_block = (Cell *)malloc(MAX_BLOCK_SIZE * sizeof(Cell));

	int num_cells_in_secondary = 0;
	int num_cells_in_main = 0;

	long file_loc = 0;
	fseek(file, 0, SEEK_END);
	long file_end = ftell(file);
	rewind(file);

	while (!feof(file)) {
		// Load a block of cells from the file
		num_cells_in_main = 0;
		
		//Limit block size and make sure the reading of file was successful.
		while (num_cells_in_main < MAX_BLOCK_SIZE && fscanf(file, "%f %f %f",  
					&main_block[num_cells_in_main].x, 
					&main_block[num_cells_in_main].y, 
					&main_block[num_cells_in_main].z) == 3) {
			num_cells_in_main++;
		}
		file_loc = ftell(file); //Tracks where next main block will start.


		// Process distances within the main block
		process_block(main_block, num_cells_in_main);

		//Reads the remaining cells into a secondary block and compares to main block.
		//If multiple blocks are required to reach the end the secondary block is overwritten
		//after each iteration.
		if(file_loc < file_end){
			while(!feof(file)){
				num_cells_in_secondary = 0;
				while (num_cells_in_secondary < MAX_BLOCK_SIZE && fscanf(file,
							"%f %f %f", 
							&secondary_block[num_cells_in_secondary].x, 
							&secondary_block[num_cells_in_secondary].y, 
							&secondary_block[num_cells_in_secondary].z) == 3){
					num_cells_in_secondary++;
				}
				//Process distances within secondary block.
				process_inter_block(main_block, num_cells_in_main, secondary_block,
						num_cells_in_secondary);

			}

			fseek(file, file_loc, SEEK_SET);

		}
	}

	// Output the distances and their counts
	output_distances(distance_count);

	// Free memory
	free(main_block);
	free(secondary_block);
	fclose(file);
	return 0;
}

