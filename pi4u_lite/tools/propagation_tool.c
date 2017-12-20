/*
 *  propagation_tool.c
 *  Pi4U
 *
 *  Created by Panagiotis Hadjidoukas on 29/6/15.
 *  Copyright 2015 ETH Zurich. All rights reserved.
 *
 */

#include "run_ssm.c" 
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>
#include <torc.h>
#include <math.h>

void taskfun( int *seed, double *res, int *info)
{
	

	printf("Executing task (%d)\n", info[0]);

	*res = run_ssm( seed, (void *)NULL, info);
	
	return;
}



int main(int argc, char *argv[])
{
	int i;

	torc_register_task(taskfun);

	torc_init(argc, argv, MODE_MS);

	srand(time(NULL));

	if( argc<2 )
		printf("\nUsage: simulate_all 100\n");


	int N = atoi(argv[1]);
	double res[N];

	for (i = 0; i < N; i++) {
		int info[1], seed[2];
		
		info[0] = i;
		seed[0]	= rand();
		seed[1]	= rand();
		
		torc_create(-1, taskfun, 3,
			2, MPI_INT, CALL_BY_COP,
			1, MPI_DOUBLE, CALL_BY_RES,
			1, MPI_INT, CALL_BY_COP,
			seed, &res[i], info);
	}

	torc_waitall();


	/*
	for (i = 0; i < t; i++) {
		printf("RESULT %03d: %10.4f %10.4f %10.4e %10.4e %10.4lf\n", i, TP[i][0], TP[i][1], ref[i], res[i], fabs(ref[i]-res[i]));
	}
	*/


	torc_finalize();

	return 0;
}

