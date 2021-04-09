/***************************************************************************
 *
 *   File        : main.c
 *   Student Id  : 763078
 *   Name        : Mustafa Awni
 *
 ***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include "tasks.h"

int main(int argc, char *argv[]) {
	
	/* TODO: Parse Command Line Arguments
	DONOT explicitly set arguments to filenames */
	char* q2_file = NULL;
	char* q4_file = NULL;
	char* q5_file = NULL;
	double xo;
	char* q6_file = NULL;
	
	q2_file = argv[1];
	q4_file = argv[2];
	q5_file = argv[3];
	xo = (double) (float) atoi(argv[4]);
		//convert argument to int. cast to float. cast to double
	//printf("xo = %f\n", xo);
	q6_file = argv[5];
	
	//valgrind ./exec in_shock.csv in_linalsys.csv in_interp.csv 4 in_waveeqn.csv
	
	/* TODO: Add timing for each task and output running time in ms */
	// timing code adapted from main.c in Ass1
	struct timeval start;
	struct timeval stop;
    
	/* Question 2 */
	gettimeofday(&start, NULL);
	shockwave(q2_file);
	gettimeofday(&stop, NULL);
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("QUESTION 2:  %.2f milliseconds\n", elapsed_ms);
	
	
	/* Question 4 */
	gettimeofday(&start, NULL);
	linalgbsys(q4_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("QUESTION 4:  %.2f milliseconds\n", elapsed_ms);
	
	/* Question 5 */
	gettimeofday(&start, NULL);
	interp(q5_file,xo);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("QUESTION 5:  %.2f milliseconds\n", elapsed_ms);
	
	/* Question 6 */
	gettimeofday(&start, NULL);
	waveeqn(q6_file);
	gettimeofday(&stop, NULL);
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0;
	elapsed_ms += (stop.tv_usec - start.tv_usec) / 1000.0;
	printf("QUESTION 6:  %.2f milliseconds\n", elapsed_ms);
    
	return (EXIT_SUCCESS);
}
