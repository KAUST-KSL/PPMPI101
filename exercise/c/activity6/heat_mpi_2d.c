/* 
 * 2D decomposition of the domain a 2D grid representing a
 * metal plate.
 * N : 			Number of points such that NxN square grid will be created as a domain
 * MAX_ITER: 	Maximum number of iterations to run
 * TOL:			Minimum delta T such that below this value no change is deemed and the simulation needs to stop
 * MAX_TEMP:	Temperature on the boundary
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define N			27
#define MAX_ITER	4000
#define TOL			1e-4
#define MAX_TEMP	100.0

float** allocate(float**);
void initialize(float**);
void halo_update(float**);
void write_grid(int, float**);



// general validity check
int is_grid_decomposible (int nprocs){
	// condition 1 : N should be exactly divisible by sqrt(nprocs)
	// condition 2 : nprocs should be a perfect square i.e. sqrt(nprocs) is a whole number
	bool cond_1=true , cond_2=true;


	int procs_1D = (int) (sqrt((double) nprocs));
	// Check first condition
	if ( N % procs_1D != 0 ){
			cond_1 = false;
	}
	// Check second condition
	if ( nprocs % procs_1D != 0 ){
			cond_2 = false;
	}
	if (( cond_1 == true ) && ( cond_2 == true)){
		return 0;
	}
	else {
		return 1;
	}
}

// Global rank and number of MPI processes
int grank, gprocs;
int globalN, localN;
// Local rank and number of MPI processes
int lrank, lprocs;
// Local rows and columns along with ghost region
int rows, cols;
float **T_old, **T_new;
// variables to store virtual topology information
int ndims = 2;			// Since we are decomposing grid in 2D
int dims[2] = { 0, 0 };  	// placeholder for dimensions
int mycoords[2] = { 0, 0 };	// placeholder for Cartesian coordinates
int period[2] = { 0, 0 };	// placeholder for periodicity information
int reorder = 0;

int x_axis = 1, y_axis = 0;
int right, left, up, down;	//placeholder for ranks of nearest neighbours

MPI_Comm cartcomm;		// Declaration of new communicator

int main(int argc, char *argv[]) {

	int i, j;


	// Initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &grank);
	MPI_Comm_size(MPI_COMM_WORLD, &gprocs);

	if ( is_grid_decomposible(gprocs) != 0){
		if (grank == 0){
			fprintf(stderr,"Grid decomposition will fail.\n");
			fprintf(stderr,"nprocs = %d  	---	Number of MPI Processes\n",gprocs);
			fprintf(stderr,"N      = %d  	---	Number of grid points in each dimension\n",N);
			fprintf(stderr,"nprocs should be a perfect square (e.g. 1,4,9,16...)\n"
					"Also, N should be exactly divisible by sqrt(nprocs)\n");
			MPI_Abort(MPI_COMM_WORLD,911);
		}
	}
	/* It may be a good idea to create a Cartesian Communicator
	 * Here are some hints:
	 * 1.	Get a good guess on dimensions
	 * 2.	Create a Cartesian communicator and call it cartcomm
	 * 3.	Query rank and its coordinates in new communicator
	 * 4.	Query who are the neighbours
	 */


	/* Domain decomposition: decide work for each rank by figuring out localN
	 * hint : globalN/number of MPI ranks
	 */

	/* Added ghost region for each MPI rank */
	rows = localN + 2;
	cols = localN + 2;

	T_old = allocate(T_old);
	T_new = allocate(T_new);

	initialize(T_old);
	initialize(T_new);
	int iter = 0;
	float dT = MAX_TEMP, dT_local = MAX_TEMP;

	while (dT > TOL && iter <= MAX_ITER) {

		halo_update(T_old);

		// Evaluate temperature in inner domain
		for (i = 1; i < rows - 1; i++) {
			for (j = 1; j < cols - 1; j++) {
				T_new[i][j] = 0.25
						* (T_old[i - 1][j] + T_old[i + 1][j] + T_old[i][j - 1]
								+ T_old[i][j + 1]);
			}
		}

		// Figure out the maximum local error
		dT = 0.0;
		dT_local = 0.0;
		for (i = 1; i < rows - 1; i++) {
			for (j = 1; j < cols - 1; j++) {
				dT_local = fmaxf(fabsf(T_new[i][j] - T_old[i][j]), dT_local);
				T_old[i][j] = T_new[i][j];
			}
		}

		/* Communicate to everyone the maximum of local errors from all ranks
		 * hint: use MPI_Allreduce with  MPI_MAX operation
		 */

		iter++;
	}


	write_grid(iter,T_old);


	// Printing results
	if (lrank == 0) {
		if ((iter - 1) == MAX_ITER)
			printf("Reached maximum iterations %d. Error = %2.4f\n", iter - 1,
					dT);
		else
			printf("Converged in %d iterations with and error of %2.4f\n",
					iter - 1, dT);
	}

	MPI_Finalize();
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

	if (mycoords[0] == dims[0] - 1) {
		for (j = 0; j < cols; j++)
			T[rows - 1][j] = MAX_TEMP;
	}
	if (mycoords[1] == 0) {
		for (i = 0; i < rows; i++)
			T[i][0] = MAX_TEMP;
	}

}

void halo_update(float **T) {

	/*
	 * Send and receive ghost regions to the 4 neighbors
	 * Hint:
	 * 1. Create temporary arrays for send and receive buffers for each side
	 * 2. load buffers
	 * 3. communicate synchronous or asynchronous, its your preference
	 * 4. unload data in receive buffer to the respective rows and columns of 2D array T.
	 */
}


void write_grid(int iter, float **T) {
	int i=0, j=0,indx=0;

	int root = 0;
	//Gather information about decomposition from all process onto Root process
	int total_elements = localN * localN;
	int* recv_numelems;
	int* disp;
	if (lrank == root) {
		recv_numelems = (int*) malloc(lprocs * sizeof(int));
		disp = (int*) malloc(lprocs * sizeof(int));
	}
	// tell Writer how many elements each rank is going to send
	MPI_Gather(&total_elements, 1, MPI_INT, recv_numelems, 1, MPI_INT, root,cartcomm);

	// disp is an array of offsets or displacements
	// i.e. start index where to receive vector from each rank (useful only on master)
	if (lrank == root) {
		disp[0] = 0;
		for (i = 1; i < lprocs; i++){
			disp[i] = disp[i-1]+recv_numelems[i];
		}
	}
	// Now prepare for gathering the actual data
	float *sendbuf = (float*) malloc(total_elements * sizeof(float));
	if (sendbuf == NULL){
		MPI_Abort(MPI_COMM_WORLD,911);
	}

	for (i = 1; i < rows - 1; i++){
		for (j = 1; j < cols - 1; j++){
			sendbuf[indx] = T[i][j];
			indx++;
		}
    }
    float* recv_buffer;
    if (lrank == root){
    	recv_buffer = (float*) malloc(globalN*globalN * sizeof(float));
    	if (recv_buffer == NULL)
    		MPI_Abort(MPI_COMM_WORLD,911);
    }
	MPI_Gatherv(sendbuf, total_elements, MPI_FLOAT, recv_buffer, recv_numelems,disp, MPI_FLOAT, root, cartcomm);

	if (lrank == root) {
		char fname[256];
		sprintf(fname, "grid_t%d.txt", iter);
		FILE *fw = fopen(fname, "w+");
		indx=0;

		int x,y;
		indx=0;
		for (y=0; y < (globalN*globalN); y+=(dims[1]*localN*localN)){
				for(i=0; i<localN*localN; i+=localN){
					for(x =0; x < dims[0]*localN*localN; x+=localN*localN){
						for(j=0;j<localN;j++){
							indx=y+i+x+j;
							fprintf(fw,"%2.5f ",recv_buffer[indx]);
						}
					}
					fprintf(fw,"\n");
				}
			}
		fclose(fw);

	}

}
