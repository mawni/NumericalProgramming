/***************************************************************************
 *
 *   File        : tasks.c
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
#include <assert.h>

/* programming is fun */
/* welcome to my submission!! Thanks for reading */

/* apologies in advance for abundance of comments. programming is fun but confusing */

/* to compile and test at home:
	gcc -Wall -std=c99 *.c -o exec -lm
	./exec in_shock.csv in_linalsys.csv in_interp.csv 4 in_waveeqn.csv
	
for submission:
login to dimefox with 
	ssh mawni@dimefox.eng.unimelb.edu.au
go to ass 2 directoy:
	cd ENGR30003/A2
check files w/: 
	ls
remove current tasks.c file when you're in the "cd ENGR30003/A2" directory:
	rm tasks.c
logout of dimefox:
	logout
copy new tasks.c file over to dimefox:
	scp tasks.c mawni@dimefox.eng.unimelb.edu.au:~/ENGR30003/A2
then proceed with submission as instructed in Ass2 pdf
	gcc -Wall -std=c99 *.c -o exec -lm
	valgrind ./exec in_shock.csv in_linalsys.csv in_interp.csv 4 in_waveeqn.csv
then submit
	submit ENGR30003 A2 *.c *.h *.pdf
	verify ENGR30003 A2 > feedback.txt
	nano feedback.txt
		note you have to press ctrl+X to close the txt
*/

/* Constants */
#define PI 3.14159265358979323846
#define THETA_MAX_M5 45
	//M=5.0, from Q2.3a this was approximated to be less than 45, so I choose 45 as the max 
	//	(I know in reality it will be lower)
#define THETA_ZERO_BETA_U 90
	//found in Q2.1
#define ACCURACY_EPSILON 0.0000001
    //epsilon for accuracy. I've chosen this num because outputs are to 6dp
//#define DELTAT 0.05
	// delta_t for Q6 timestepping
	// FALSE!! We get delta t from CFL
#define MAXTIME 0.2
	//for Q6, it's the last timestep to be analysed. 0<=t<=0.2

/* Q2 shockwave struct */
typedef struct {
	double M;	//mach number
	double theta; //angle in degrees (not radians)
		//I STORED AS DOUBLE, BUT REMEMBER IT MUST BE OUTPUTTED AS INTEGER
	double beta_l;
	double beta_u;
		//lower and upper beta. degrees (not radians)
	double gamma; //ratio of specific heats
} shockwave_t;

/* Q4 vector struct */
typedef struct{
	double a;
	double b;
	double c;
	double Q;
		//these are taken directly from A2 pdf
	double x;
		//basically x_i for vector {X}
} vector_t;

/* Q5 spline struct */
	//technically the last data point doesn't have a spline, but no worries
typedef struct{
	double x;
	double f_x;
		//this is just f(x)
	double a,b,c,d;
		//the constants/coefficients for each spline
		//S_i(x) = a_i + b_i*(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
	double tri_diag_a, tri_diag_b, tri_diag_c, tri_diag_Q;
		//these are for a tridiagonal system that will be solved to get all c_i's for the splines
		//completely unrelated to the double a,b,c,d's defined in this struct
		//I used tri_diag_a|b|c|Q to be consistent with Q4
} spline_t;

/* Q6 wave point struct */
typedef struct{
	double x;
	double fx;
		//f(x)
} wavept_t;



/* taken from ENGR3K3>L2>filo_io.c */
FILE* safe_fopen(const char* path, const char* mode)
{
    FILE* fp = fopen(path, mode);
    if (fp == NULL) {
        perror("file open error.");
        exit(EXIT_FAILURE);
    }
    return fp;
}

/* taken from workshops */
/*
void* safe_malloc(size_t bytes)
{
    void* ptr = malloc(bytes);
    if (ptr == NULL) {
        printf("error mallocing %lu bytes\n", bytes);
        exit(EXIT_FAILURE);
    }
    return ptr;
}*/

/* Newton raphson eqns Q2.3 */
/* Remember that c uses radians */
double f_shockwave(double M, double theta, double beta, double gamma) {
	return (2*(cos(beta)/sin(beta))*(pow(M,2)*pow(sin(beta),2)-1)/(pow(M,2)*(gamma+cos(2*beta))+2)-tan(theta));
}

double fprime_shockwave(double M, double beta, double gamma) {
	//derivative of f_shockwave. I calculated by hand
	return ((-2*pow(1/sin(beta),2)*(pow(M,2)*pow(sin(beta),2)-1)+pow(M,2)*sin(2*beta)*2*(cos(beta)/sin(beta)))*(pow(M,2)*(gamma+cos(2*beta))+2)+2*pow(M,2)*sin(2*beta)*2*(cos(beta)/sin(beta))*(pow(M,2)*pow(sin(beta),2)-1))/(pow((pow(M,2)*(gamma+cos(2*beta))+2),2));
}

double newtonraph(double beta_i, double M, double theta, double gamma) {
	return beta_i-f_shockwave(M, theta, beta_i, gamma)/fprime_shockwave(M, beta_i, gamma);
}


void shockwave(const char* q2_file)
{
    //printf("shockwave() - IMPLEMENT ME!\n");
    
    /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 0. define struct has M, theta, beta_l, beta_u, gamma
    //		### not necessary but i'll make it anyway
    //			I'll write beta_l and beta_u to file as they are calculated 
    // 1. open in_shock.csv file, 
    //		1.1. read first 5 lines and do nothing
    //		1.2. read until read failure for how many M values there are. count the number of lines
    // 2. close file. reset to line 1
    // 3. malloc 1D array of M values to be analysed
    // 4. reopen in_shock.csv file
    //		4.1. take in data from first few lines then input all M values into 1D array
    // 5. close file
    // 6. Q2.3a do newton raphson method for M=5, theta=20
    //		6.1. ####WHAT TO DO ABOUT INITIAL GUESSES??####
    //		6.2. print final beta_l and beta_u values. check with ass2 pdf
    // 7. malloc array of structs
    //		7.1. i'll need to define a struct that contains theta, beta_l, beta_u
    // 8. Q2.3b do newton raphson method for M=5, 0<=theta<=theta_max
    //		8.1. loop through theta = 0 until theta = theta_max(??) or until theta = 90 deg
    //			8.1.1. do newton raphson for theta = 0 using gormula from Q2.1 in Ass2
    // 9. maybe print these numbers to a .csv but i don't want that csv to be submitted
    //		print it and then take the values into excel and make the graph. put into excel
    // 10. free malloc array of structs for (7)
    // 11. create out_shock.csv for outputting part data
    //		11.1. print initial default line 1 "M,theta,beta_lower,beta_upper\n"
    // 12. Q2.3c do newton-raphson loop over all M values in array of Ms (3)
    //		12.1. for each M, run through theta = 0 until theta = theta_max and calculate beta_l and beta_u
    //			12.1.1. fprintf to csv ("%f,%d,%f,%f\n",M,theta,beta_l,beta_u) 
    // 13. close out_shock.csv
    // free malloc array for M (3)
    
    
    //### Open in_shock.csv file, 
    FILE *f;
	f = safe_fopen(q2_file,"r");
    
	
    //### read first 5 lines and do nothing
    	//there's a lot of debug code commented out. Had issues reading file but got it eventally
    fscanf(f,"M,theta,beta_l,beta_u,gamma\n");
    fscanf(f,"%*f,%*f,%*f,%*f,%*f\n");
    /*  
	double temp;
    fscanf(f,"%lf,",&temp);
    printf("%f\n",temp);
    fscanf(f,"%lf,",&temp);
    printf("%f\n",temp);
    fscanf(f,"%lf,",&temp);
    printf("%f\n",temp);
    fscanf(f,"%lf,",&temp);
    printf("%f\n",temp);
    fscanf(f,"%lf",&temp);
    printf("%f\n",temp);
    */
    //double temp;
    //char temp_c;
    	//debug stuff
    fscanf(f,"M\n");
    //if (fscanf(f,"%c",&temp_c)!=1){printf("success at line 145\ntemp_c = [%c]\n",temp_c);};
    	//1 data item is matched and stored stored, so returned value will be 1
    	//fscanf(f,"M"); wasn't working for me. It would never successfully read the M, so I did %c instead
    //if (fscanf(f,"5.0")!=0){printf("error at line 148\n");};
    //fscanf(f,"M");
    fscanf(f,"%*f\n");
    fscanf(f,"M");
    //fscanf(f,"%lf\n",&temp);
    //printf("temp = %f\n", temp);
    /* %*f skips the variable */
    
    //### read until read failure for how many M values there are. count the number of lines
    int M_val_ctr = 0;
    while (fscanf(f,"\n%*f") == 0){
    	M_val_ctr += 1;
    }
    //printf("For Q2.3c, there are %d M values to analyse\n", M_val_ctr);
    	//should come out as 6. SUCCESS!
    
    //### close file
    fclose(f);	
    
    //### malloc 1D array of M values to be analysed
	double* M_vals = NULL;
	M_vals = (double*) malloc((M_val_ctr)*sizeof(double));
		// c stores array from M_vals[0] until M_vals[M_val_ctr-1]
	assert(M_vals != NULL);
    
    //### Re-open file
	f = safe_fopen(q2_file,"r");
	if (fseek(f, 0, SEEK_SET)) {
		//reset file so I can read again from beginning
		perror("file reset error");
    	exit(EXIT_FAILURE); 
	}
    
	//### Defining values for Q2.3a, Q2.3b
	double M_a, M_b;
	double theta_a; //degrees
	double beta_l_a, beta_u_a;
		//lower and upper beta. degrees (not radians). these are initial guesses
	double gamma; //ratio of specific heats
	
    //### take in data from first few lines then input all M values into 1D array
    fscanf(f,"M,theta,beta_l,beta_u,gamma\n");
    fscanf(f,"%lf,%lf,%lf,%lf,%lf\n", &M_a, &theta_a, &beta_l_a, &beta_u_a, &gamma);
    fscanf(f,"M\n");
    fscanf(f,"%lf\n", &M_b);
    fscanf(f,"M");
    for (int i = 0; i < M_val_ctr; i++){
    	// loop through until expected M_val_ctr M values are read
    	fscanf(f,"\n%lf", &M_vals[i]);
    		//store M number into the array
    	//printf("M number #%d = %lf\n", i+1, M_vals[i]);
    }
    
    //### close file
    fclose(f);
    
    beta_l_a = beta_l_a*(PI/180);
    beta_u_a = beta_u_a*(PI/180);
    theta_a = theta_a*(PI/180);
    	//180 degrees not really a magic number in my opinion. Pretty obvious
    	//convert to radians	
    
    //double function_test = f_shockwave(5.0, 0.0, PI/2, 1.4);
    //printf("For M = 5, theta = 0, beta_u = 90 deg, gamma = 1.4 --> f(beta) = %f\n", function_test);
    	//checking my function is correct. And it is. 
    //double derivative_test = fprime_shockwave(5.0, PI/2, 1.4);
    //printf("For M = 5, beta_u = pi/2 radians, gamma = 1.4 --> f'(beta) = %f\n", derivative_test);
    	//checking my derivative is typed up correct. And it is.
    	//https://www.wolframalpha.com/input/?i=derivative+of+%282*%28cos%28x%29%2Fsin%28x%29%29*%285%5E2*sin%28x%29%5E2-1%29%2F%285%5E2*%281.4%2Bcos%282*x%29%29%2B2%29-tan%280%29%29+where+x+%3D+pi%2F2+radians
	
    //### Q2.3a do newton raphson method for M=5, theta=20
    int maxiter = 100;
    	//i'll stop iterating when 100 loops of newtonraph() have been calculated
    //printf("Q2.3a. For M = %f, theta = %f rads, gamma = %f\n", M_a, theta_a, gamma);
    //printf("Iteration 0, beta_l = %f, f(beta_l) = %f \n", beta_l_a, f_shockwave(M_a, theta_a, beta_l_a, gamma));
    for (int i=1; i<(maxiter+1); i++){
    	beta_l_a = newtonraph(beta_l_a, M_a, theta_a, gamma);
    		//iterate newtonraph beta value
    	//printf("Iteration %d, beta_l = %f, f(beta_l) = %f \n", i, beta_l_a, f_shockwave(M_a, theta_a, beta_l_a, gamma));
    	if (fabs(f_shockwave(M_a, theta_a, beta_l_a, gamma)) < ACCURACY_EPSILON){
    		//we've found an accurate enough rot
    		i = maxiter+1;
    	}
    }
    
    //printf("--\nIteration 0, beta_u = %f, f(beta_l) = %f \n", beta_u_a, f_shockwave(M_a, theta_a, beta_u_a, gamma));
    for (int i=1; i<(maxiter+1); i++){
    	beta_u_a = newtonraph(beta_u_a, M_a, theta_a, gamma);
    	//printf("Iteration %d, beta_u = %f, f(beta_u) = %f \n", i, beta_u_a, f_shockwave(M_a, theta_a, beta_u_a, gamma));
    	if (fabs(f_shockwave(M_a, theta_a, beta_u_a, gamma)) < ACCURACY_EPSILON){
    		i = maxiter+1;
    	}
    }
    
    

    //######## Q2.3b do newton raphson method for M=5, 0<=theta<=theta_max
    /* it's been commented out because it generates a csv file for my own analysis. 
    //To just comment out that file printing section would be a pain ._.
    
    
    //int theta_b;
    	//will be analysed from 0<=theta<=theta_max 
    double beta_l_b, beta_u_b;
    	//beta values for part (b)
    beta_l_b = asin(1/M_b);
		//Eqn from Q2.1 for when theta = 0. Stored in radians. It will be my initial guess for theta = 1 deg
	beta_u_b = THETA_ZERO_BETA_U*(PI/180);
		//from Q2.1 for when theta = 0
	
	FILE* fp = safe_fopen("q2_3b.csv", "w");
		// 'q2_3b.csv' doesnt exist yet, so we create it with write mode 'w'
	    //### maybe print these numbers to a .csv but i don't want that csv to be submitted. 
	    //I'll make the csv once and use it to make the graph for Q2.3b
    if (fprintf(fp, "theta,beta_l,beta_u\n") < 0){
		perror("file write error");
        exit(EXIT_FAILURE);  
	}
    fprintf(fp, "0.0,%f,%f\n", beta_l_b*(180/PI), beta_u_b*(180/PI));
    	//print the data for when theta = 0. we already know the solutions
	
    int beta_u_solve_chk = 0; //0 if unsolved, 1 if solved
    int beta_l_solve_chk = 0;
    for (int theta_b=1; theta_b<(THETA_MAX_M5+1); theta_b++){
    	beta_u_solve_chk = 0;
    	beta_l_solve_chk = 0;
    	for (int i=1; i<(maxiter+1); i++){
    		if (beta_u_solve_chk==0){
    			beta_u_b = newtonraph(beta_u_b, M_b, theta_b*(PI/180), gamma);
    				//previous solution for beta is the initial guess for the next one
    				//done because of the trends seen in Q2.2a and Fig 2 in A2 pdf
    				if(fabs(f_shockwave(M_b, theta_b*(PI/180), beta_u_b, gamma)) < ACCURACY_EPSILON){
    					//if root found accurately, stop iterating
    					beta_u_solve_chk = 1;
    				};
    		};
    		if (beta_l_solve_chk==0){
    			beta_l_b = newtonraph(beta_l_b, M_b, theta_b*(PI/180), gamma);
    			if(fabs(f_shockwave(M_b, theta_b*(PI/180), beta_l_b, gamma)) < ACCURACY_EPSILON){
    					beta_l_solve_chk = 1;
    			};
    		};
    		
    		if (beta_l_solve_chk == 1 && beta_u_solve_chk == 1){
    			i = maxiter+1;
    			fprintf(fp, "%d,%f,%f\n", theta_b, beta_l_b*(180/PI), beta_u_b*(180/PI));
    			//### print it and then take the values into excel and make the graph. put into excel
    		};
    	};
    };
    	
    //following section not necessary/
    // 7. malloc array of structs
    //		7.1. i'll need to define a struct that contains theta, beta_l, beta_u
    //		8.1. loop through theta = 0 until theta = theta_max(??) or until theta = 90 deg?
    //free malloc array of structs for (7)
    
    fclose(fp);
    
    */

    
    //############## Q2.3c
    
    //### create out_shock.csv for outputting part data
	FILE* fp_c = safe_fopen("out_shock.csv", "w");
		// 'out_shock.csv' doesnt exist yet, so we create it with write mode 'w'
	
	//### print initial default line 1 "M,theta,beta_lower,beta_upper\n"
    if (fprintf(fp_c, "M,theta,beta_lower,beta_upper\n") < 0){
		perror("file write error");
        exit(EXIT_FAILURE);  
	};
	
	double beta_l_c, beta_u_c;
    	//beta values for part (b)
    int beta_u_solve_chk; //0 if unsolved, 1 if solved
	int beta_l_solve_chk;
	
    //### for each M, run through theta = 0 until theta = 90 and calculate beta_l and beta_u
    // 	my algorithm works such that if no solutions exist for a theta, nothing is outputted
    // 	this means that by just picking theta = 90, I'm playing it quite safe
    //	thinking intuitively, it makes no sense for theta_max to even be close to 90, so it should be fine
	for (int k=0; k < M_val_ctr; k++){
		beta_l_c = asin(1/M_vals[k]);
			//Eqn from Q2.1 for when theta = 0. Stored in radians. 
			//It's my answer for when theta=0, and my initial guess for theta = 1 deg
		beta_u_c = THETA_ZERO_BETA_U*(PI/180);
			//theta_zero_beta_u = 90
			//from Q2.1 for when theta = 0
		fprintf(fp_c, "%f,0,%f,%f\n", M_vals[k], beta_l_c*(180/PI), beta_u_c*(180/PI));
    		//print the data for when theta = 0 deg. we already know the solutions
    		
		beta_u_solve_chk = 0; //0 if unsolved, 1 if solved
		beta_l_solve_chk = 0;
		
		for (int theta_c=1; theta_c<(THETA_ZERO_BETA_U+1); theta_c++){
			//to make it as dynamic as possible, I've set the algorithm to run until theta = 90
			//this is always the maximum, so it's a safe move. explained in another comment around line 351
			beta_u_solve_chk = 0;
			beta_l_solve_chk = 0;
			for (int i=1; i<(maxiter+1); i++){
				if (beta_u_solve_chk==0){
					//if soln not yet found, try again
					beta_u_c = newtonraph(beta_u_c, M_vals[k], theta_c*(PI/180), gamma);
					//previous solution for beta is the initial guess for the next one
					//done because of the trends seen in Q2.2a and Fig 2 in A2 pdf
					if(fabs(f_shockwave(M_vals[k], theta_c*(PI/180), beta_u_c, gamma)) < ACCURACY_EPSILON){
						//if root found accurately, stop iterating
						beta_u_solve_chk = 1;
					};
				};
				if (beta_l_solve_chk==0){
					//if soln not yet found, try again
					beta_l_c = newtonraph(beta_l_c, M_vals[k], theta_c*(PI/180), gamma);
					if(fabs(f_shockwave(M_vals[k], theta_c*(PI/180), beta_l_c, gamma)) < ACCURACY_EPSILON){
						beta_l_solve_chk = 1;
					};
				};
				
				if ((beta_l_solve_chk == 1) && (beta_u_solve_chk == 1) && (beta_l_c >= 0)){
					//if both solutions found, print to csv. (beta_l_c >= 0) is an added precaution that solns are realisable
					//fprintf to csv required data
					fprintf(fp_c, "%f,%d,%f,%f\n", M_vals[k], theta_c, beta_l_c*(180/PI), beta_u_c*(180/PI));
					//printf("M=%f,theta=%d,beta_l=%f,beta_u=%f\n",M_vals[k], theta_c, beta_l_c*(180/PI), beta_u_c*(180/PI));
					i = maxiter+1;
				};
			};
			//printf("theta = %d, beta_l soln found? = %d, beta_u soln found? = %d\n", theta_c, beta_l_solve_chk, beta_u_solve_chk);
			if((beta_l_solve_chk == 0) || (beta_u_solve_chk == 0) || (beta_l_c < 0)){
				//if by whatever incremented theta_c value, if no soln found, means we've reached theta_max
				//the (beta_l_c<0) is an added precaution so that if by 100 iterations, 
					//value becomes negative, it's not realisable anymore, 
					//and all subsequent solutions are detached
				//stop finding solutions and end loop of incrementing theta_c
				theta_c = THETA_ZERO_BETA_U+1;
				//printf("no solutions found. theta = %d\n", theta_c);
			};
		};
    };
    
    fclose(fp_c);
    free(M_vals);
    
}

void linalgbsys(const char* q4_file)
{
    //printf("linalgbsys() - IMPLEMENT ME!\n");
    
    /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 0. define a vector struct has a,b,c,Q
    // 1. open q4_file AKA in_linalsys.csv file, 
    //		1.1. read first line and do nothing
    //		1.2. read lines until read failure for how many vectors there are. count the number of lines
    // 2. close file
    // 3. malloc 1D array of structs to be analysed
    // 4. reopen in_linalsys.csv file. reset to line 1 of file
    //		4.1. input all vectors into 1D array until all lines read
    // 5. close file
    // 6. Q4 do the back substritution method
    //		6.1. if i start calculations at last array element (row N), then can move back procedurally
    // 7. create out_linalsys.csv for outputting solution vector {X} (6dp)
    //		7.1. loop through array of structs and print each x_i
    // 8. close out_linalsys.csv
    // 9. free malloc array for M (3)
    

    //### open q4_file AKA in_linalsys.csv file 
    FILE *f;
	f = safe_fopen(q4_file,"r");
	
    fscanf(f,"a,b,c,q");
    	//read first line and do nothing
    	
    //### read lines until read failure for how many vectors there are. count the number of lines
    int ctr_N = 0;
    while (fscanf(f,"\n%*f,%*f,%*f,%*f") == 0){
    	ctr_N += 1;
    }
    //printf("For Q4, there are %d rows of vector to read in\n", ctr);
    	//should come out as 6. SUCCESS!
    
    //### close file
    fclose(f);
    
    // 3. malloc 1D array of structs to be analysed
    vector_t* vectors = NULL;
    vectors = (vector_t*) malloc(ctr_N*sizeof(vector_t));
		// c stores array from vectors[0] until vectors[ctr_N-1]
	assert(vectors != NULL);
    
    //### reopen in_linalsys.csv file. reset to line 1 of file
    f = safe_fopen(q4_file,"r");
    if (fseek(f, 0, SEEK_SET)) {
		//reset file so I can read again from beginning
		perror("file reset error");
    	exit(EXIT_FAILURE); 
	}
    
    //### input all vectors into array until all lines read
    fscanf(f,"a,b,c,q");
    for (int i = 0; i < ctr_N; i++){
    	// loop through until expected ctr_N rows are read
    	fscanf(f,"\n%lf,%lf,%lf,%lf", &vectors[i].a, &vectors[i].b, &vectors[i].c, &vectors[i].Q);
    		//store numbers into the array
    	//printf("Vector row %d --> %lf,%lf,%lf,%lf\n", i+1, vectors[i].a, vectors[i].b, vectors[i].c, vectors[i].Q);
    };
    
    //### close file
    fclose(f);
    
    //### starting from row 2, calculate all a_i and Q_i
    	//for row 1, a_i and Q_i remain unchanged
    for (int i=1 ; i<ctr_N ; i++){
    	vectors[i].a = vectors[i].a - vectors[i].c * (vectors[i-1].b / vectors[i-1].a);
    	vectors[i].Q = vectors[i].Q - vectors[i].c * (vectors[i-1].Q / vectors[i-1].a);
    		//eqns taken from A2 pdf
    };
    
    //### Q4 do the back substritution method
    vectors[ctr_N-1].x = vectors[ctr_N-1].Q / vectors[ctr_N-1].a; 
    
    for (int i = ctr_N-2; i > -1; i--){
    	//starting from second last element of vectors. back substitution
    	vectors[i].x = (vectors[i].Q - vectors[i].b * vectors[i+1].x) / vectors[i].a; 
    		//using formula given in A2 pdf
    };
    
    //### create out_linalsys.csv for outputting solution vector {X} (6dp)
    FILE* f_c = safe_fopen("out_linalsys.csv", "w");
		// 'out_linalsys.csv' doesnt exist yet, so we create it with write mode 'w'
	
	//### print initial default line"
    //ass 2 pdf doesn't specify what file format should be
    if (fprintf(f_c, "x\n") < 0){
		perror("file write error");
        exit(EXIT_FAILURE);  
	};
    
	fprintf(f_c,"%lf", vectors[0].x);
		//print first x_i, i=0
    //### loop through array of structs and print each x_i
    for (int i = 1; i < ctr_N; i++){
    	fprintf(f_c,"\n%lf", vectors[i].x);
    	//printf("Vector row %d --> %lf,%lf,%lf,%lf\n", i+1, vectors[i].a, vectors[i].b, vectors[i].c, vectors[i].Q);
    }
    
    //### close out_linalsys.csv
    fclose(f_c);
    
    //### free malloc array
    free(vectors);
    
}

void interp(const char* q5_file, const double xo)
{
    //printf("interp() - IMPLEMENT ME!\n");
    
    /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 0. define a struct? x_i and f(x_i)
    //		0.1. cubic spline S_i(x) = a_i + b_i*(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
    //		0.2. also store the spline info. a,b,c,d
    //		0.3. note that the last x_i shouldn't have a corresponding spline(?)
    //			0.3.1. ####OR DO I NEED A PHANTOM SPLINE
    // 1. open q5_file AKA in_interp.csv file, 
    //		1.1. read first line and do nothing
    //		1.2. read lines until read failure for how many given data points there are. count the number of lines
    // 2. close file
    // 3. malloc 1D array of structs to be analysed?
    // 4. reopen in_interp.csv file. reset to line 1 of file
    //		4.1. input all data poitns into 1D array until all lines read
    // 5. close file
    // 6. Q5 do the interpolation
    // 7. create out_interp.csv for outputting 'xo' and 'f(xo)'
    //		7.1. print first line "xo,f(xo)\n"
    //		7.2. calculate f(xo) solution(s) for xo
    //			7.3. if there are 2, need to output both solns
    // 8. close out_interp.csv
    // 9. free malloc array for M (3)
    
    
    
    //### open q5_file AKA in_interp.csv file,
    FILE *f;
    f = safe_fopen(q5_file,"r");
    
    fscanf(f, "x,f(x)\n");
    	//read first line and do nothing
    int pt_ctr = 0;
    while (fscanf(f,"%*f,%*f") == 0){
    	//read lines until read failure for how many given data points there are. count the number of lines
    	pt_ctr += 1;
    }
    //printf("For Q5, there are %d data points given\n", pt_ctr);
    	//should come out as 52. SUCCESS!
    
    //### close file
    fclose(f);
    
    //### malloc 1D array of structs to be analysed
    spline_t* splines = NULL;
    splines = (spline_t*) malloc(pt_ctr*sizeof(spline_t));
    assert(splines!=NULL);
    
    //### reopen in_interp.csv file. reset to line 1 of file
    f = safe_fopen(q5_file,"r");
    if (fseek(f, 0, SEEK_SET)) {
		//reset file so I can read again from beginning
		perror("file reset error");
    	exit(EXIT_FAILURE); 
	}
	
	fscanf(f, "x,f(x)");
    //### input all data poitns into 1D array until all lines read
    for (int i = 0; i < pt_ctr ; i++){
    	fscanf(f,"\n%lf,%lf", &splines[i].x, &splines[i].f_x);
    		//store numbers into the array
    	//printf("Data point %d --> %.15lf,%.15lf\n", i, splines[i].x, splines[i].f_x);
    		//everything is read correctly for our given file!
    };
    
    //### close file
    fclose(f);
    
    //########## DO THE INTERPOLATION ##########
    
    //### set all a_i = f(x_i)
    for (int i = 0; i < pt_ctr ; i++){
    	splines[i].a = splines[i].f_x;
    };
    
    
    //### natural spline because we know no other info. Set endpoint second derivatives to 0
    	// c_0=c_n=0
    splines[0].c = 0;
    splines[pt_ctr-1].c = 0;
    	// Since i=0 to n --> pt_ctr = n+1 --> n = pt_ctr-1
    
    //### from lecture 18, a tridiagonal system exists for relationship btwn a_i's and c_i's
    	// it looks like [A]{X]={C} 
    		//[A] is a tridiagonal matrix with relationships btwn the x_i values
    		//{X} is a vector for all the c_i's
    		//{C} is a vector containing relationships btwn x_i's and a_i's
    
    //### first we calculate all the non-zero elements in [A]. This is like matrix form (6) from A2 pdf
    	//more specifically though, it's from Lecture 18, Slide 69 (L18,S69). That shows the matrix format
    	//in spline_t, the matrix elements are defined as tri_diag_a & tri_diag_b & tri_diag_c
    	//i'm using the same format used in Q4 of this project
    
    splines[0].tri_diag_a = 1;
    splines[0].tri_diag_b = 0;
    	//first row of matrix according to L18,S69 is as above
    splines[pt_ctr-1].tri_diag_a = 1;
    splines[pt_ctr-1].tri_diag_c = 0;
    	//last row of matrix according to L18,S69 is as above
    for (int i=1; i<pt_ctr-1; i++){
    	//not calculating first and last row
    	splines[i].tri_diag_c = splines[i].x - splines[i-1].x;
    	splines[i].tri_diag_a = 2*(splines[i+1].x- - splines[i-1].x);
    	splines[i].tri_diag_b = splines[i+1].x - splines[i].x;
    		/* We look at the matrix in L18,S69 and consider the top row as
    		row 0, and last row as row n (like how n is the last data point for 
    		interpolation. Then we compare that matrix with (6) from Q4 in the 
    		Ass2 pdf. By linking these two together, we can associate a pattern
    		for each tri_diag_c|a|b and create these eqns I've got above. Note
    		that I took the eqns from L18,S69 and simplified the h_i into x_i's */
    }
    //Now we've successfully generated the matrix in L18,S69
    
    
    //### now calculate all elements in {C}. This is just like {Q} in Q4
    	// i'll do these calculations using the eqns from L18,S70
    splines[0].tri_diag_Q = 0;
    splines[pt_ctr-1].tri_diag_Q = 0;
    	//first and last element are as above according to L18,S70
    for (int i=1; i<pt_ctr-1; i++){
    	//not calculating first and last row
    	splines[i].tri_diag_Q = (3/(splines[i+1].x - splines[i].x))*(splines[i+1].a-splines[i].a) + (3/(splines[i].x - splines[i-1].x))*(splines[i-1].a-splines[i].a);
    };
    
    
    //Now all parts of tridiagonal system are ready to do
    //### Rewrite the matrix into the form (7) from Q4 in Ass2 pdf
        	// note for row 0, tri_diag_a and tri_diag_Q remain unchanged
    for (int i=1 ; i<pt_ctr ; i++){
    	splines[i].tri_diag_a = splines[i].tri_diag_a - splines[i].tri_diag_c * (splines[i-1].tri_diag_b / splines[i-1].tri_diag_a);
    	splines[i].tri_diag_Q = splines[i].tri_diag_Q - splines[i].tri_diag_c * (splines[i-1].tri_diag_Q / splines[i-1].tri_diag_a);
    		//eqns taken from A2 pdf
    };	
    //Now we're ready to complete the back substitution
    
    
    //### do the back substitution to solve for {X} (the vector of c_i's)
    splines[pt_ctr-1].c = splines[pt_ctr-1].tri_diag_Q / splines[pt_ctr-1].tri_diag_a; 
    
    
    for (int i = pt_ctr-2; i > -1; i--){
    	//starting from second last element of vectors. back substitution
    	splines[i].c = (splines[i].tri_diag_Q - splines[i].tri_diag_b * splines[i+1].c) / splines[i].tri_diag_a; 
    		//using formula given in A2 pdf
    };
    //Now by this point, we have all a_i's and all c_i's for the splines
    
    //printf("currently on line 721\n");
    
    //### calculate all b_i's for the splines
    for (int i = 0; i<pt_ctr-1; i++){
    	//last point has no spline
    	//printf("line 720. i = %d\n", i);
    	splines[i].b = (splines[i+1].a - splines[i].a)/(splines[i+1].x - splines[i].x) - ((splines[i+1].x - splines[i].x)*(2*splines[i].c + splines[i+1].c))/3;
    		//eqn taken from L18,S67. Eqn (55)
    }
    
    //printf("currently on line 729\n");
    
    //### calculate all d_i's for the splines
    for (int i = 0; i<pt_ctr-1; i++){
    	//last point has no spline
    	//printf("line 729. i = %d\n", i);
    	splines[i].d = (splines[i+1].c - splines[i].c)/(3*(splines[i+1].x - splines[i].x));
    }
    
    //FINALLY all splines have been calculated
    //printf("currently on line 737\n");
    
    //### create out_interp.csv for outputting 'xo' and 'f(xo)'
    FILE* f_c = safe_fopen("out_interp.csv", "w");
    if (fprintf(f_c, "xo,f(xo)") < 0){
		perror("file write error");
        exit(EXIT_FAILURE);  
	};
    
	double f_xo;
		//the interpolated value for xo using our spline
    //### take 'xo' and see whenever it's found within an interval, output
    for (int i=0; i<(pt_ctr-1); i++){
    	if ((splines[i].x <= xo) && (xo <= splines[i+1].x)){
    		//x lies in spline i
    		f_xo = splines[i].a + splines[i].b*(xo - splines[i].x) + splines[i].c*pow((xo - splines[i].x),2) + splines[i].d*pow((xo - splines[i].x),3);
    		fprintf(f_c, "\n%f,%f", xo, f_xo);
    	} else if ((splines[i+1].x <= xo) && (xo <= splines[i].x)){
    		//since the function loops around, i'm also testing the interval where x_i > x_(i+1), spline moves right to left
    		f_xo = splines[i].a + splines[i].b*(xo - splines[i].x) + splines[i].c*pow((xo - splines[i].x),2) + splines[i].d*pow((xo - splines[i].x),3);
    		fprintf(f_c, "\n%f,%f", xo, f_xo);
    	}
    	//printf("Spline %d: =%f + %f*(4 - $$4) + %f*(4 - $$4)^2 + %f*(4 - $$4)^3\n",i, splines[i].a, splines[i].b, splines[i].c, splines[i].d);
    		//for excel
    	//printf("Spline %d: %f + %f*(x - x_i) + %f*(x - x_i)^2 + %f*(x - x_i)^3\n",i, splines[i].a, splines[i].b, splines[i].c, splines[i].d);
    		//for testing
    	//S_i(x) = a_i + b_i*(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3
    }
    
    
    //### close out_interp.csv
    fclose(f_c);
    
    
    //### free malloc array (3)
    free(splines);
    
    //this problem should've been worth more than 4 marks!! -_-
    
}

void waveeqn(const char* q6_file)
{
    //printf("heateqn() - IMPLEMENT ME!\n");
    
    /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 0. have a struct for wavept_t
    // 1. open q6_file AKA in_interp.csv file, 
    //		1.1. read first line "c,Nx,CFL,out_iter\n"
    //		1.2. read second line and take in data
    // 2. close file
    // 3. calculate number of points N_x+1 according to delta_x=1/N_x over 0<=x<=1
    // 4. malloc 3 1D array of x_i's & func values. f_i_n, f_i_n+1, f_i_n+0.5 (intermediate) where i=0,1,2,...,N_x
    // 5. Q6 solve the differential eqn
    //		5.1. set initial conditions for f(x,t=0)
    // 6. create out_waveeqn_1U.csv & out_waveeqn_2C.csv
    //		6.1. fprintf outout
    // 7. close files from (6)
    // 8. free malloc array for M (4)
    
    //### open q6_file AKA in_interp.csv file,
    FILE *f;
    f = safe_fopen(q6_file,"r");
    fscanf(f, "c,Nx,CFL,out_iter\n");
    
    double c, CFL;
    int Nx, out_iter;
    //### read second line and take in data
    fscanf(f, "%lf,%d,%lf,%d", &c, &Nx, &CFL, &out_iter);
    //printf("c = %lf, Nx = %d, CFL = %lf, out_iter = %d\n", c, Nx, CFL, out_iter);
    //### close file
    fclose(f);
    
    //### calculate deltax = 1.0/Nx over 0<=x<=1
    double deltax = 1.0/Nx;
    //printf("delta x = %lf\n", deltax);
    	//comes out as deltax = 0.005. correct for our given Nx
    
    //malloc 3 1D arrays of x_i's & func values. f_i_n, f_i_n+1, f_i_n+0.5 (intermediate) where i=0,1,2,...,N_x
    wavept_t* wavepts_fn = NULL;
    wavept_t* wavepts_fn1 = NULL;
    wavept_t* wavepts_fn05 = NULL;
    	//fn = f_i_n, fn1 = f_i_n+1, fn05 = f_i_n+0.5
    wavepts_fn = (wavept_t*) malloc((Nx+1)*sizeof(wavept_t));
    wavepts_fn1 = (wavept_t*) malloc((Nx+1)*sizeof(wavept_t));
    wavepts_fn05 = (wavept_t*) malloc((Nx+1)*sizeof(wavept_t));
    	//note that in total, there are Nx+1 points
    assert(wavepts_fn!=NULL);
    assert(wavepts_fn1!=NULL);
    assert(wavepts_fn05!=NULL);
    
    //### initialise all x_i values according to the deltax and Nx
    for (int i = 0; i<(Nx+1); i++){
    	wavepts_fn[i].x = i*deltax;
    	//printf("x_%d = %f\n", i, wavepts_fn[i].x);
    	wavepts_fn1[i].x = i*deltax;
    	wavepts_fn05[i].x = i*deltax;
    }
               
    //FILE* fn = safe_fopen("out_fn_initialcond.csv","w");
    	//file for printing intitial conditions for f(x,t=0)
    //fprintf(fn, "t=0\n");
    //fprintf(fn, "x,fx");
    	
    //### set initial conditions for f(x,t=0) according to A2 pdf Q6
    for (int i = 0; i<(Nx+1); i++){
    	if ((0 <= wavepts_fn[i].x) && (wavepts_fn[i].x < 0.125)){
    		//0.125 isn't really a magic number, it's just given A2 pdf
    		wavepts_fn[i].fx = 0;
    	} else if ((0.125 <= wavepts_fn[i].x) && (wavepts_fn[i].x <= 0.375)){
    		//again, 0.375 isn't a magic number. it's a very obvious one actually =P
    		wavepts_fn[i].fx = 0.5*(1-cos(8*PI*(wavepts_fn[i].x-0.125)));
    			//eqn given in A2 pdf Q6
    	} else{
    		// 0.375<x<=1
    		wavepts_fn[i].fx = 0;
    	}
    	//fprintf(fn, "\n%lf,%lf", wavepts_fn[i].x, wavepts_fn[i].fx);
    	//printf("%lf,%lf\n", wavepts_fn[i].x, wavepts_fn[i].fx);
    }
    //fclose(fn);
    
    //### calculate delta_t
    double deltat = (CFL*deltax)/c;
    //printf("delta t = %f, CFL = %f\n", deltat, CFL);
    int timesteps = (MAXTIME / deltat)+1;
    //printf("total number of timesteps = %d\n", timesteps);
    	//it's an int because it only makes sense for the delta t to equally split up the time interval
    	//a non-integer amount of timesteps should not happen
    //printf("amount of timesteps starting from and including zero = %d\n", timesteps);
    	//should equal 5. success!
    
    	
    FILE* fnt = safe_fopen("out_waveeqn_1U.csv","w");
		//for writing out desired timestep solution
	fprintf(fnt, "x,f(x)");	
	
	double RHS_fn, RHS_fn05;	// for storing RHS(fn) and RHS(fn+0.5)
	
    //### make a loop to go over each timestep and within loop, runge-kutta 2nd order will run
	for (int i = 0; i<(timesteps-1); i++){
		//using timestep 0 to get next timestep. Means the last timestep is calculated using the previous timestep 
		//current time = i*deltat
		
		//#######################Q6 run runge-kutta 2nd order and solve the differential eqn
		
		/* ~~~~~~~~~~~~~~~~~~ */
		/* FIRST ORDER SCHEME */
		/* ~~~~~~~~~~~~~~~~~~ */
		
		//#### find fn05 from fn
		//loop through all x points
		for(int j = 0; j<(Nx+1); j++){
			//run for x points in interval for fn
			//calculate RHS's for fn
			if (j==0){
				//use boundary stencil when j=0
				RHS_fn = (wavepts_fn[1].fx - wavepts_fn[0].fx)/deltax;
			} else {
				//use upwind scheme for every next j=1,2,...,Nx
				RHS_fn = (wavepts_fn[j].fx-wavepts_fn[j-1].fx)/deltax;
			}	
			//use eqn (11) from A2 pdf to get fn05
			wavepts_fn05[j].fx = wavepts_fn[j].fx + deltat*RHS_fn;
			
		}
					
		//#### find fn1 from fn05
		for(int j = 0; j<(Nx+1); j++){
			//run for x points in interval for fn+0.5
			
			//### calculate RHS's for fn
			if (j==0){
				//use boundary stencil when j=0
				RHS_fn = (wavepts_fn[1].fx - wavepts_fn[0].fx)/deltax;
			} else {
				//use upwind scheme for every next j=1,2,...,Nx
				RHS_fn = (wavepts_fn[j].fx-wavepts_fn[j-1].fx)/deltax;
			}
			
			//### calculate RHS's for fn+0.5
			if (j==0){
				//use boundary stencil when j=0
				RHS_fn05 = (wavepts_fn05[1].fx - wavepts_fn05[0].fx)/deltax;
			} else {
				//use upwind scheme for every next j=1,2,...,Nx
				RHS_fn05 = (wavepts_fn05[j].fx-wavepts_fn05[j-1].fx)/deltax;
			}	
			//use eqn (11) from A2 pdf to get fn+1
			wavepts_fn1[j].fx = wavepts_fn[j].fx + (deltat/2)*(RHS_fn-RHS_fn05);
			
		}		
		//### NOW FN+1 HAS BEEN SOLVED. WE'VE FOUND THE NEXT TIMESTEP SOLN
		//### print out desired timestep solution
		
		for (int j = 0; j<(Nx+1); j++){
			if ((out_iter==0)&&(i==0)){
				//if we want to print out starting time and we're currently looking at it, print current fn
				fprintf(fnt, "\n%f,%f", wavepts_fn[j].x, wavepts_fn[j].fx);
			} else if (/*i==out_iter*/i*deltat==0.1){
				//else we want to print whatever other timestep out_iter equals
				fprintf(fnt, "\n%f,%f", wavepts_fn1[j].x, wavepts_fn1[j].fx);
			}
		}
		//#### reset. set wavepts_fn[i] = wavepts_fn1[i]
			//run for all x
			//do this so the next loop looks at the solved wavepts to get new wavepts_fn1
		for (int j = 0; j<(Nx+1); j++){
			wavepts_fn[j].fx = wavepts_fn1[j].fx;
		}	
		
	}
    fclose(fnt);	
    
    
    /* ~~~~~~~~~~~~~~~~~~~ */	
	/* SECOND ORDER SCHEME */
	/* ~~~~~~~~~~~~~~~~~~~ */
	FILE* fnt2C = safe_fopen("out_waveeqn_2C.csv","w");
		//for writing out desired timestep solution
	fprintf(fnt2C, "x,f(x)");	
	
	//double RHS_fn, RHS_fn05;	// for storing RHS(fn) and RHS(fn+0.5)
	//were already declared above
	
    //### make a loop to go over each timestep and within loop, runge-kutta 2nd order will run
	for (int i = 0; i<(timesteps-1); i++){
		//using timestep 0 to get next timestep. Means the last timestep is calculated using the previous timestep 
		//current time = i*deltat
		//printf("i = %d, current time i*deltat = %f\n",i, i*deltat);
		//#######################Q6 run runge-kutta 2nd order and solve the differential eqn
		
		//#### find fn05 from fn
		//loop through all x points
		for(int j = 0; j<(Nx+1); j++){
			//run for x points in interval for fn
			//calculate RHS's for fn
			if (j==0){
				//use boundary stencil when j=0
				RHS_fn = (wavepts_fn[1].fx - wavepts_fn[0].fx)/deltax;
			} else if (i==Nx){ 
				//use boundary stencil when j=0
				RHS_fn = (wavepts_fn[Nx].fx - wavepts_fn[Nx-1].fx)/deltax;
			} else {
				//use central scheme for every next i=1,2,...,Nx-1
				RHS_fn = (wavepts_fn[j+1].fx-wavepts_fn[j-1].fx)/(2*deltax);
			}	
			//use eqn (11) from A2 pdf to get fn05
			wavepts_fn05[j].fx = wavepts_fn[j].fx + deltat*RHS_fn;	
		}
		//printf("RHS_fn = %f\n", RHS_fn);
		//#### find fn1 from fn05
		for(int j = 0; j<(Nx+1); j++){
			//run for x points in interval for fn+0.5
			
			//### calculate RHS's for fn
			if (j==0){
				//use boundary stencil when j=0
				RHS_fn = (wavepts_fn[1].fx - wavepts_fn[0].fx)/deltax;
			} else if (j==Nx){
				//use boundary stencil when j=Nx
				RHS_fn = (wavepts_fn[Nx].fx - wavepts_fn[Nx-1].fx)/deltax;
			} else {
				//use central scheme for every next i=1,2,...,Nx-1
				RHS_fn = (wavepts_fn[j+1].fx-wavepts_fn[j-1].fx)/(2*deltax);
			}
			
			//### calculate RHS's for fn+0.5
			if (j==0){
				//use boundary stencil when j=0
				RHS_fn05 = (wavepts_fn05[1].fx - wavepts_fn05[0].fx)/deltax;
			} else if(j==Nx){
				//use boundary stencil when j=Nx
				RHS_fn05 = (wavepts_fn05[Nx].fx - wavepts_fn05[Nx-1].fx)/deltax;
			} else {
				//use central scheme for every next i=1,2,...,Nx-1
				RHS_fn05 = (wavepts_fn05[j+1].fx-wavepts_fn05[j-1].fx)/(2*deltax);
			}	
			//use eqn (11) from A2 pdf to get fn+1
			wavepts_fn1[j].fx = wavepts_fn[j].fx + (deltat/2)*(RHS_fn-RHS_fn05);
			
		}		
		
		//### NOW FN+1 HAS BEEN SOLVED. WE'VE FOUND THE NEXT TIMESTEP SOLN
		//### print out desired timestep solution
		
		for (int j = 0; j<(Nx+1); j++){
			if ((out_iter==0)&&(i==0)){
				//if we want to print out starting time and we're currently looking at it, print current fn
				fprintf(fnt2C, "\n%f,%f", wavepts_fn[j].x, wavepts_fn[j].fx);
			} else if (/*i==out_iter*/(i*deltat)==0.1){
				//else we want to print whatever other timestep out_iter equals
				fprintf(fnt2C, "\n%f,%f", wavepts_fn1[j].x, wavepts_fn1[j].fx);
			}
		}
		//#### reset. set wavepts_fn[i] = wavepts_fn1[i]
			//run for all x
			//do this so the next loop looks at the solved wavepts to get new wavepts_fn1
		for (int j = 0; j<(Nx+1); j++){
			wavepts_fn[j].fx = wavepts_fn1[j].fx;
		}	
		
	}
    fclose(fnt2C);
		
    
    
    //### free malloc arrays
    free(wavepts_fn);
    free(wavepts_fn1);
    free(wavepts_fn05);
    
    
}
