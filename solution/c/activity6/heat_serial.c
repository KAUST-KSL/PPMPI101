/*
 * 1D decomposition of the domain in the column dimension
 * such that each MPI task will get an equal number of rows.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N			27
#define MAX_ITER	4000
#define TOL			1e-4
#define MAX_TEMP	100.0

float** allocate(float**);
void initialize(float**);
void write_grid(int, float**);

int rows, cols;
float **T_old, **T_new;

int main(int argc, char *argv[]) {

	int i, j;

	rows = N + 2, cols = N + 2;

	T_old = allocate(T_old);
	T_new = allocate(T_new);

	initialize(T_old);
	initialize(T_new);

	int iter = 0;
	float dT = MAX_TEMP;

	while (dT > TOL && iter <= MAX_ITER) {
		for (i = 1; i < rows - 1; i++) {
			for (j = 1; j < cols - 1; j++) {
				T_new[i][j] = 0.25
						* (T_old[i - 1][j] + T_old[i + 1][j] + T_old[i][j - 1]
								+ T_old[i][j + 1]);
			}
		}
		dT = 0.0;
		for (i = 1; i < rows - 1; i++) {
			for (j = 1; j < cols - 1; j++) {
				dT = fmaxf(fabsf(T_new[i][j] - T_old[i][j]), dT);
				T_old[i][j] = T_new[i][j];
			}
		}
		iter++;
	}

	write_grid(iter, T_old);

	if ((iter - 1) == MAX_ITER)
		printf("Reached maximum iterations %d. Error = %2.4f\n", iter, dT);
	else
		printf(
				"Converged in %d iterations with and error of %2.4f\n", iter, dT);

	return 0;
}

float** allocate(float **T) {
	int i = 0, j = 0;

	// Allocate including some extra for ghost regions

	T = (float**) malloc(rows * sizeof(float*));
	if (T != NULL) {
		T[0] = (float*) malloc(rows * cols * sizeof(float));
	}
	if (T[0] != NULL) {
		for (i = 0; i < rows; i++)
			T[i] = (*T + i * cols);
	}
	return T;
}
void initialize(float **T) {
	int i, j;

	for (i = 0; i < rows; i++)
		for (j = 0; j < cols; j++)
			T[i][j] = 0.75 * MAX_TEMP;

	for (j = 0; j < cols; j++)
		T[rows - 1][j] = MAX_TEMP;

	for (i = 0; i < rows; i++)
		T[i][0] = MAX_TEMP;

}

void write_grid(int iter, float **T) {
	int i = 0, j = 0, indx = 0;

	int root = 0;
	//Gather information about decomposition from all process onto Root process
	int total_elements = N * N;

	char fname[256];
	sprintf(fname, "output_serial_t%d.txt", iter);
	FILE *fw = fopen(fname, "w+");

	for (i = 1; i < rows - 1; i++) {
		for (j = 1; j < cols - 1; j++) {
			fprintf(fw, "%2.5f ", T[i][j]);
		}
		fprintf(fw, "\n");
	}

	fclose(fw);

}
