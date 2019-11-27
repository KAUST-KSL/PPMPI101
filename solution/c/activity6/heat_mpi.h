int is_grid_decomposible (int nprocs){
	// condition 1 : N should be exactly divisible by sqrt(nprocs)
	// condition 2 : nprocs should be a perfect square i.e. sqrt(nprocs) is a whole number
	int cond_1=false , cond_2=false;


	int procs_1D = (int) (sqrt((double) nprocs);
	// Check first condition
	if ( N % procs_1D != 0 )
			cond_1 = True;
	// Check second condition
	if ( sqrt((double) nprocs) % procs_1D != 0 )
			cond_2 = True;

	if (( cond_1 == True ) && ( cond_2 == True )){
		return 0;
	}
	else {
		return 1;
	}
}
