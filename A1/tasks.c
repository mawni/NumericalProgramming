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
	gcc -Wall -std=c99 *.c -o flow -lm
	./flow flow_data.csv 10
	
for submission:
login to dimefox with 
	ssh mawni@dimefox.eng.unimelb.edu.au
go to ass 1 directoy:
	cd ENGR30003/A1
check files w/: 
	ls
remove current tasks.c file when you're in the "cd ENGR30003/A1" directory:
	rm tasks.c
logout of dimefox:
	logout
copy new tasks.c file over to dimefox:
	scp tasks.c mawni@dimefox.eng.unimelb.edu.au:~/ENGR30003/A1
then proceed with submission as instructed in Ass1 pdf
	gcc -Wall -std=c99 *.c -o flow -lm
	valgrind ./flow flow_data.csv 10
then submit
	submit ENGR30003 A1 *.c *.h
	verify ENGR30003 A1 > feedback.txt
	nano feedback.txt
		note you have to press ctrl+X to close the txt
*/

/* Constants */
//#define FILENAME "./flow_data.csv"
	//needs to be read in
/* Total amount of points checked by opening flow_data.csv in excel*/
#define TOTAL_FLOW_PTS 419629
	// Didn't see anywhere on the assignment that says this isn't allowed!
//#define GRID_RESOLUTION 10
	//needs to be read in 
#define MIN_X -15
#define MAX_X 85
#define MIN_Y -20
#define MAX_Y 20
	/* min and max values given in Ass 1 pdf */
	
// from workshop 3	
#define BST_SUCCESS 1
#define BST_FAILURE 0

#define BST_PREORDER 0
#define BST_INORDER 1
#define BST_POSTORDER 2


/* Task 1 flow pt struct */
typedef struct {
	double rho;
	double u;
	double v;
	double x;
	double y;
} flow_t1;

/* Task 2 flow pt struct */
typedef struct {
	double rho;
	double u;
	double v;
	double x;
	double y;
	double score;
	int k;
		/* Counts how many flow points are included in the coarse grid cell. */
		/* Used for calculating averages */
} coarse_flow_t2;

/* Task 3 flow pt struct */
typedef struct {
	double rho;
	double u;
	double v;
	double x;
	double y;
	double rho_u;
		//rho*u
} flow_t3;

/* Task 4 flow pt struct */
typedef struct {
	double rho;
	double u;
	double v;
	double x;
	double y;
	double w;
		//vorticity
} flow_t4;


/* adapted from workshop 2 */
/* shoutout to Chitrarth to his amazing solutions to the workshops */
/* node type */
typedef struct node node_t;
	//recursion. node within nodes. node_t declared before used
struct node {
	void* data;
    //flow_t3 data;
    node_t* next;
    //linear search will be done in task 3, so don't need double linked list
};
/* linked list type */
typedef struct {
    int num_elements;
    node_t* head;
    node_t* tail;
    void (*del)(void*);
} list_t;



/* from workshop 3 */
/* bst node type */
typedef struct node_bst node_bst_t;

struct node_bst {
    void* data;
    node_bst_t* left;
    node_bst_t* right;
};
/* bst type */
typedef struct {
    int num_elements;
    node_bst_t* root;
    void (*del)(void*);
    int (*cmp)(const void*, const void*);
} bst_t;



/* from workshop 2 */
/* remove node at the front of the list */
void* list_pop_front(list_t* list)
{
    assert(list != NULL);
    assert(list->num_elements > 0);
    node_t* old;
    assert(list->head != NULL);
    old = list->head;
    list->head = list->head->next;
    void* d = old->data;
    free(old);
    list->num_elements--;
    if (list->num_elements == 0) {
        list->head = NULL;
        list->tail = NULL;
    }
    return d;
}

/* from workshop 2 */
/* add node add the front of the list */
void list_push_front(list_t* list, void* d)
{
    assert(list != NULL);
    node_t* new = (node_t*)malloc(sizeof(node_t));
    assert(new);
    new->data = d;
    new->next = list->head;
    	//make new node's 'next' pointer point to whatever node the current head points to
    	//new has been pushed into the front
    list->head = new;
    	//make current head point to 'new' node
    if (list->tail == NULL)
        list->tail = new;
    		// if tail points to NULL, make it point to the 'new' node.
    		// ###DONT ACTUALLY UNDERSTAND WHY###
    list->num_elements++;
}


/* from workshop 2 */
/* create a new linked list structure */
list_t* list_new(void (*delfunc)(void*))
	//list_new() returns a pointer to list_t
{
    list_t* list;
    list = (list_t*)malloc(sizeof(list_t));
    	//allocate memory as normal
    assert(list != NULL);
    list->head = NULL;
    list->tail = NULL;
    	//standard procedure. head and tail nodes = NULL
    list->num_elements = 0;
    list->del = delfunc;
    return list;
}

/* from worshop 2 */
/* free all memory associated with a list */
void list_free(list_t* list)
{
    assert(list != NULL);
    while (list->num_elements) {
        void* elem = list_pop_front(list);
        list->del(elem);
    }
    free(list);
}

int doublecmp(double a, double b){
    if (a>b){
    	return 1;
    } else if (a<b){
    	return -1;
    } else {
    	return 0;
    }
}

//cmp_t3(double a, double b)
int cmp_t3(const void * a, const void * b){
	flow_t3 *flow_A = (flow_t3 *)a;
		//declare pointer *flow_A that points to a flow_t3 struct
		//cast 'a' as pointer for flow_t3
		//set ptr flow_A to cast ptr a
		//REMEMBER this function is made for ONLY qsort(). But somehow worked for bst and list
		
	flow_t3 *flow_B = (flow_t3 *)b;
    if ((flow_A->rho_u)>(flow_B->rho_u)){
    	return 1;
    	//struct element pointed by 'a' goes AFTER 'b'
    } else if ((flow_A->rho_u)<(flow_B->rho_u)){
    	return -1;
    	//struct element pointed by 'a' goes BEFORE 'b'
    } else {
    	return 0;
    }
}

/* from workshop 3 */
void no_free(void* v)
{
}

/* from workshop 3 */
/* create a new empty bst structure */
bst_t* bst_new(void (*delfunc)(void*), int (*cmpfunc)(const void*, const void*))
{
    bst_t* bst;
    bst = (bst_t*)malloc(sizeof(bst_t));
    assert(bst != NULL);
    bst->root = NULL;
    bst->num_elements = 0;
    bst->del = delfunc;
    bst->cmp = cmpfunc;
    return bst;
}

/* modified from workshop 3 */
/* insert a new element into the bst */
int bst_insert(bst_t* bst, void* d)
{
    assert(bst != NULL);
    assert(d != NULL);
    node_bst_t* parent = NULL;
    node_bst_t* tmp = bst->root;
    	//make node_bst_t 'tmp' = whatever node root points to
    while (tmp) {
    	//while tmp is not NULL
        parent = tmp;
        if (bst->cmp(tmp->data, d) > 0) { // element is smaller
            tmp = tmp->left;
        }
        else if (bst->cmp(tmp->data, d) < 0) { // element is bigger
            tmp = tmp->right;
        }
        else {
            /* ALREADY EXISTS! -> do nothing and return ERROR */
            return BST_FAILURE;
        }
    }

    /* insert as child of parent */
    node_bst_t* new_node = (node_bst_t*)malloc(sizeof(node_bst_t));
    assert(new_node != NULL);
    new_node->data = d;
    new_node->left = NULL;
    new_node->right = NULL;
    if (parent != NULL) {
        if (bst->cmp(parent->data, d) > 0) { // element is smaller
            assert(parent->left == NULL);
            parent->left = new_node;
        }
        else {
            assert(parent->right == NULL);
            parent->right = new_node;
        }
    }
    else {
        assert(bst->root == NULL);
        bst->root = new_node;
    }
    bst->num_elements++;

    return BST_SUCCESS;
}

/* modified from workshop 3 */
void perfect_insert(bst_t* bst, flow_t3* array, int low, int high)
{
    if (low <= high) {
    	// Choose root from array and insert
    	// Recursively do the same on left and right (1)
        int mid = low + (high - low) / 2;
        flow_t3* ptr = array + mid;
        	//array is a pointer to an array of flow_t3 structs
        	//mid is just the middle element number
        	//ptr points to a flow_t3, it is = array+mid because we are jumping from first element to mid
        bst_insert(bst, ptr);
        perfect_insert(bst, array, low, mid - 1);
        perfect_insert(bst, array, mid + 1, high);
    }
}


/* from workshop 3 */
/* free all memory assocated with a subtree */
void bst_free_subtree(bst_t* bst, node_bst_t* n)
{
    assert(bst != NULL);
    if (n) {
        bst_free_subtree(bst, n->left);
        bst_free_subtree(bst, n->right);
        bst->del(n->data);
        free(n);
        bst->num_elements--;
    }
}
/* free all memory associated with a bst */
void bst_free(bst_t* bst)
{
    assert(bst != NULL);
    bst_free_subtree(bst, bst->root);
    free(bst);
}

/* taken from workshop 3 */
/* dont need. didn't finish modifying it. kept here just in case
int make_unique(flow_t3* array, int n)
{
    int dest = 0;
    int itr = 1;
    while (itr != n) {
        if (array[dest] != array[itr]) {
            array[++dest] = array[itr];
        }
        itr++;
    }
    return dest+1;
}*/

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

/* function to check if 'a' perfectly divides 'b' */
/* will be used in task 2 */
// ###UPDATE### Now using 'fmodf(_,_)' instead. returns remainder of float division
/*
int perfect_division(float a, float b){
	if (a % b == 0){
		// no remainder
		return 1; 
			// 'b' is perfect divisor
    }
    else {
    	// there is remainder
    	return 0;
    		// 'b' NOT perfect divisor
    }
} 
*/

void maxfluxdiff(const char* flow_file)
{
    //printf("maxfluxdiff() - IMPLEMENT ME!\n");
    
    /* Open flow_data.csv file */
    FILE *f;
	f = safe_fopen(flow_file,"r");
	//fscanf(f,"%d",&size); 
	
	// Allocate memory 
	flow_t1* flows = NULL;
	flows = (flow_t1*) malloc((TOTAL_FLOW_PTS+1)*sizeof(flow_t1));
		// it's TOTAL_FLOW_PTS+1 because c starts arrays from a[0] but I want to start analysis from a[1]
	assert(flows != NULL);
	
    /* for (i=0;i<5;i++) --> analyses i= 0 to 4, total 5 loops */
    
    /*initial scan to disregard first line*/
    int i = 1;
    fscanf(f,"rho,u,v,x,y%lf,%lf,%lf,%lf,%lf",&flows[i].rho, &flows[i].u, &flows[i].v, &flows[i].x, &flows[i].y); 
    
    /* %*f skips the y variable */
    
    
    //printf("%f,%f,%f,%f,%f\n", flows[i].rho, flows[i].u, flows[i].v, flows[i].x, flows[i].y);
    
    /* the elements in array for desired flux values */
    int i_maxflux_u = 0, i_minflux_u = 0, i_maxflux_v = 0, i_minflux_v = 0;
    /* the actual flux values. initialising as first elements */
    double maxflux_u = flows[i].rho*flows[i].u;
    double minflux_u = flows[i].rho*flows[i].u;
    double maxflux_v = flows[i].rho*flows[i].v;
    double minflux_v = flows[i].rho*flows[i].v;
    
    /* for processing */
    double flux_u, flux_v; 
    int result_u_max, result_u_min, result_v_max, result_v_min; 
    
    /* scan in all data */
    for (i=2;i<(TOTAL_FLOW_PTS+1);i++){
    	fscanf(f,"%lf,%lf,%lf,%lf,%lf",&flows[i].rho, &flows[i].u, &flows[i].v, &flows[i].x, &flows[i].y);
    	
    	/*####### I DON'T THINK I NEED AN ARRAY OF STRUCTS. #######
    		Why not continously fscanf values and overwrite the previous ones. 
    		Just have array of 4 structs, max & min u, max & min v.
    		200 IQ if I do this. 
    	*/
    	
    	if (flows[i].x > 20){
    		flux_u = flows[i].rho*flows[i].u;
    		flux_v = flows[i].rho*flows[i].v;
    		/* do the comparisons for min and max rho*u and rho*v */
    		result_u_max = doublecmp(flux_u, maxflux_u);
    		result_u_min = doublecmp(flux_u, minflux_u);
    		result_v_max = doublecmp(flux_v, maxflux_v); 
    		result_v_min = doublecmp(flux_v, minflux_v); 
    		
    		/* set new flux values if conditions are met */
    		if (result_u_max > 0){
    			maxflux_u = flows[i].rho*flows[i].u;
    			i_maxflux_u = i;
    		} else if (result_u_min < 0){
    			minflux_u = flows[i].rho*flows[i].u;
    			i_minflux_u = i;
    		} 
    		
    		if (result_v_max > 0){
    			maxflux_v = flows[i].rho*flows[i].v;
    			i_maxflux_v = i;
    		} else if (result_v_min < 0){
    			minflux_v = flows[i].rho*flows[i].v;
    			i_minflux_v = i;
    		}
    	}
    }
    
    
    //i = 419629; /* the last element in the array */
    //printf("last rho = %f\n",flows[i-1].rho);
    /* Note that i = 20000 means row 20002 in flow_data.csv when opened in excel */
    
    fclose(f);
    
    /* Print the values. for testing */
    /*
    printf("rho,u,v,x,y\n");
    printf("%f,%f,%f,%f,%f\n", flows[i_maxflux_u].rho, flows[i_maxflux_u].u, flows[i_maxflux_u].v, flows[i_maxflux_u].x, flows[i_maxflux_u].y);
    printf("%f,%f,%f,%f,%f\n", flows[i_minflux_u].rho, flows[i_minflux_u].u, flows[i_minflux_u].v, flows[i_minflux_u].x, flows[i_minflux_u].y);
    printf("%f,%f,%f,%f,%f\n", flows[i_maxflux_v].rho, flows[i_maxflux_v].u, flows[i_maxflux_v].v, flows[i_maxflux_v].x, flows[i_maxflux_v].y);
    printf("%f,%f,%f,%f,%f\n", flows[i_minflux_v].rho, flows[i_minflux_v].u, flows[i_minflux_v].v, flows[i_minflux_v].x, flows[i_minflux_v].y);
    */
    
    	
    
    /* write data to file */
    /* modified from ENGR3K3>L2>filo_io.c */
    FILE* fp = safe_fopen("task1.csv", "w");
        /* 'task1.csv' doesnt exist yet, so we create it with write mode 'w' */
    /* remember to first include the file pointer, then begin printing stuff */
    if (fprintf(fp, "rho,u,v,x,y\n")
    	&& fprintf(fp, "%f,%f,%f,%f,%f\n", flows[i_maxflux_u].rho, flows[i_maxflux_u].u, flows[i_maxflux_u].v, flows[i_maxflux_u].x, flows[i_maxflux_u].y)
    	&& fprintf(fp, "%f,%f,%f,%f,%f\n", flows[i_minflux_u].rho, flows[i_minflux_u].u, flows[i_minflux_u].v, flows[i_minflux_u].x, flows[i_minflux_u].y)
    	&& fprintf(fp, "%f,%f,%f,%f,%f\n", flows[i_maxflux_v].rho, flows[i_maxflux_v].u, flows[i_maxflux_v].v, flows[i_maxflux_v].x, flows[i_maxflux_v].y)
    	&& fprintf(fp, "%f,%f,%f,%f,%f\n", flows[i_minflux_v].rho, flows[i_minflux_v].u, flows[i_minflux_v].v, flows[i_minflux_v].x, flows[i_minflux_v].y) < 0){
        perror("file write error");
        exit(EXIT_FAILURE);  
    }
    fclose(fp);
        
    
    free(flows);  
}

void coarsegrid(const char* flow_file, int resolution)
{
    //printf("coarsegrid() - IMPLEMENT ME!\n");
    
    /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 1. have variable that defines grid size, figure out the borders. maybe look at multiples of delta_x or delta_y until max borders of entire grid
    // 2. malloc for 2D array of structs for coarse grid
    // 		2.1. include a counter element in each struct for num of points in a coarse grid cell
    //		2.2. malloc for 1D array of structs for scanf data
    //			2.2.1. ###DON'T NEED THIS### just directly perform (3.1)-(3.4)
    // 3. loop that scanfs data from file
    // 		3.1. if point lies on the border of a cell, ignore it
    //		3.2. checks which cell if its values lie in a coarse grid cell 
    //		3.3. sum the rho,u,v,x,y vals into respective coarse grid cell --> store this in a 2D array (2)
    //		3.4. increment counter element whenever flow pt processed
    // 4. divide sums of rho,u,v,x,y by num of points 'k'
    // 5. calculate score
    // 6. sort the data in descending order of score
    // 7. print all data to task2.csv 
    // 8. malloc free() coarse grid array
    
    
    /* Open file */
    FILE *f;
	f = safe_fopen(flow_file,"r");
	
	//#####VALGRIND ERROR CHECK#####
	//printf("currently on line 267\n");
	
	/* Allocate memory for 2D array coarse grid */ 
	coarse_flow_t2** coarse_flows = NULL;
	coarse_flows = (coarse_flow_t2**) malloc((resolution+1)*sizeof(coarse_flow_t2));
	assert(coarse_flows != NULL);	
	for (int i = 1; i < (resolution+1); i++){
		coarse_flows[i] = (coarse_flow_t2*) calloc((resolution+1),sizeof(coarse_flow_t2));
			// ####### it's resolution+1 because c starts arrays from a[0] but I want to start analysis from a[1]
			// since i've doing a bunch of += and /= later in some loops, i use calloc to initialise to 0
			// note that instead of malloc((resolution+1)*sizeof(coarse_flow_t2)), calloc is written as above w/ the comma
	}
	
	//#####VALGRIND ERROR CHECK#####
	//printf("currently on line 279\n");
    
	
	/* grid cell size values coarse flow grid */
	double delta_x = (MAX_X-MIN_X)/resolution;
	//printf("delta x = %f\n", delta_x);
	double delta_y = (MAX_Y-MIN_Y)/resolution;
	//printf("delta y = %f\n", delta_y);
	
	// -15 -05 05 15 25 35 45 55 65 75 85 
	// in the case of grid resolution 10, 11 x vals must be checked
	// 11 --> resolution+1
	
	/* declare variables for fscanf, these are overwritten with each loop */
	double rho=0;
	double u=0;
	double v=0;
	double x=0;
	double y=0;
	
	/* ignore first line*/
    fscanf(f,"rho,u,v,x,y"); 
    
	
	int i=0;
	int cell_border_check=0; // 0 = not on border, 1 = cell on border
	int x_cell_num=0;
	int y_cell_num=0;
	/* scan in all data */
	for (i=1;i<(TOTAL_FLOW_PTS+1);i++){
    	fscanf(f,"%lf,%lf,%lf,%lf,%lf",&rho, &u, &v, &x, &y);

    	/* debug
    	if (i==1){
    		printf("%f,%f,%f,%f,%f\n", rho, u, v, x, y);
    	}*/
    	
    	cell_border_check = 0;
    		/* reset the cell border check */
    		
    	//if point lies on the border of a cell, ignore it
    	if (fmodf(x-MIN_X, delta_x) == 0){
    		/* difficult to explain this. I'm changing the range of -15 to 85
    		(MIN_X to MAX_X) into 0 to 100. For the case of 
    		resolution = 10, when |-15| is	added to them, the values 
    		-15, -05, 05, ..., 85 change. They become 0, 10, ..., 90, 100. 
    		These are multiples of delta_x.	delta_x*0,1,...,resolution.
    		##WHEN THE REMAINDER OF (x-MIN_X)/delta_x = 0, THE FLOW POINT WILL
    		BE LYING ON THE BORDER OF A COARSE GRID CELL##
    		*/
    		cell_border_check = 1;
    	} else if (fmodf(y-MIN_Y, delta_y) == 0){
    		cell_border_check = 1;
    	}
    	/* Not gonna lie, pretty proud of this. */
    	
    	/* abandoned approach to check if on border. much more tedious
    	for(i=0;i<(resolution+1);i++){ 
    		// for worst case, loop will run resolution+1 (11) times
    		if (x == (MIN_X+i*delta_x)){
    			i = resolution + 1;
    				// want loop to end. 
    		} else if(y == (MIN_Y+i*delta_y)){
    			i = resolution + 1;
    		}	
    	}*/
    	
    	if (cell_border_check == 0){
    		
    		// find which coarse grid cell the flow pt is in
    		// runs for:
    		//              i = 0     1     2     3     4     5     6     7     8     9   NOT 10
    		//MIN_X+i*delta_x = -15   -5    5     15    25    35    45    55    65    75  NOT 85
    		//           cell =    1     2     3     4     5     6     7     8     9      10               
    		// don't need to look at x=85 because it's used for the cell between 75 and 85
    		
    		for(int i2=0;i2<(resolution);i2++){ 
    			if ((MIN_X+i2*delta_x) < x && x < (MIN_X+(i2+1)*delta_x)){
    				x_cell_num = i2+1; 
    			}	
    		}
    		for(int i3=0;i3<(resolution);i3++){ 
    			if ((MIN_Y+i3*delta_y) < y && y < (MIN_Y+(i3+1)*delta_y)){
    				y_cell_num = i3+1;
    				
    				/*debug. after i=418571, --> segmentation fault. 
    				i=418571 is the last valid point though, all later ones lie
    				on a coarse grid cell boundary
    				if (i>=418560){
    					printf("i = %d, x cell num = %d, y cell num = %d\n", i, x_cell_num, y_cell_num);
    				}
    				*/ 
    			}	
    		}
    		
    		/*
    		if (i>=418560){
    			printf("i = %d was assigned a cell\n", i);
    			printf("sum of u currently = %f\n", coarse_flows[x_cell_num][y_cell_num].u);
    		}*/
    		
    		/* debug
    		if(i==1060){
    			// i know from the flow_data.csv file that:
    			// point 1060 isn't on a border
    			// point 1059 isn't on a border
    			// point 1058 IS on a border
    			printf("for flow point i = %d, u = %f, x = %f, y = %f, and cell_border_check = %d\n", i, u, x, y, cell_border_check);
    			printf("x cell = %d, y cell = %d\n", x_cell_num, y_cell_num);
    		}*/
    		
    		//	sum the rho,u,v,x,y vals into respective coarse grid cell --> store this in a 2D array (2)
    		//printf("i = %d, x = %f, ", i, x);
    		coarse_flows[x_cell_num][y_cell_num].rho += rho;
    		//printf("u = %f\n", u);
    		coarse_flows[x_cell_num][y_cell_num].u += u;
    			//before error was fixed, when this line for 'u' was uncommented
    			//there was a segmentation fault! --> problem fixed by increasing the malloc by a bit
    		coarse_flows[x_cell_num][y_cell_num].v += v;
    		coarse_flows[x_cell_num][y_cell_num].x += x;
    		coarse_flows[x_cell_num][y_cell_num].y += y;
    		coarse_flows[x_cell_num][y_cell_num].k += 1;
    			//	increment counter element whenever flow pt processed into cell
    		
    	}
    	
    	/* debug
    	if(i==9){
    		// i know from the flow_data.csv file that point 9 is on border
    		printf("for flow point i = %d, x = %f, y = %f, and cell_border_check = %d\n", i, x, y, cell_border_check);
    	}
    	if(i==1){
    		// i know from the flow_data.csv file that point 1 is on border
    		printf("for flow point i = %d, x = %f, y = %f, and cell_border_check = %d\n", i, x, y, cell_border_check);
    	}*/
    	
    
    }
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 419\n");
	
    fclose(f);
    //printf("scanf loop completed. i = %d\n", i);
    	//if all items are successfully scanned, i will be printed as i=419630 aka TOTAL_FLOW_PTS+1
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 426\n");
		//many many valgrind errors after this line
		//#### FIXED!!!!!
	
    // go through all elements in 2D coarse grid array
    for (int i=1; i < (resolution+1); i++){
    	for (int j=1; j < (resolution+1); j++){
    		//divide sums of rho,u,v,x,y by num of points 'k'
    		//find averages
    		
    		//#####VALGRIND ERROR CHECK#####
    		//printf("currently on line 436, [i][j] = [%d][%d]\n",i,j );
    		
    		coarse_flows[i][j].rho /= coarse_flows[i][j].k;
    		coarse_flows[i][j].u /= coarse_flows[i][j].k;
    		coarse_flows[i][j].v /= coarse_flows[i][j].k;
    		coarse_flows[i][j].x /= coarse_flows[i][j].k;
    		coarse_flows[i][j].y /= coarse_flows[i][j].k;
    		//printf("element [i,j] aka [%d,%d] averages processed successfull\n", i,j);
    		
    		//#####VALGRIND ERROR CHECK#####
    		//printf("currently on line 446, [i][j] = [%d][%d]\n",i,j );
    		
    		//calculate score
    		coarse_flows[i][j].score = 100*sqrt((pow(coarse_flows[i][j].u,2) + 
    			pow(coarse_flows[i][j].v,2))/(pow(coarse_flows[i][j].x,2) + 
    				pow(coarse_flows[i][j].y,2)));
    	}
    }
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 449\n");
		// many valgrind errors before this check
		//####FIXED!!!!!! I was using += and /= throughout the code. Used calloc to initialise everything to 0
	
    // sort the data in descending order of score
    int k=0, l=0;
    	// [i][j] <-->
    	// [k][l] <--> [row][column]
    coarse_flow_t2* temp_coarse_flow = (coarse_flow_t2*) malloc(sizeof(coarse_flow_t2));
    //coarse_flow_t2 temp_coarse_flow;
    assert(temp_coarse_flow != NULL);
    //sorting a 2D array. Took an eon to write
	for (int j = 1; j < (resolution+1); j++) {
		//printf("j = %d\n", j);
		for (i = 1; i < (resolution+1); i++) {
			//i choose to sort in a row. [1,j]--descending-->[1,j]
			//then change row, and continue. [2,j]--descending-->[2,j]. etc
			*temp_coarse_flow = coarse_flows[i][j];
				//set the intial temp flow point
			k = i + 1;
				//for incoming initial loop, skip current element
			for (l = j; l < (resolution+1); l++) {
				//loop through y-axis
				while (k < (resolution+1)) {
					//loop through x-axis while locked on a row (y) 	
					if (coarse_flows[k][l].score > temp_coarse_flow->score) {
						// if (current element is > previous temp one)
						*temp_coarse_flow = coarse_flows[k][l];
						coarse_flows[k][l] = coarse_flows[i][j];
						coarse_flows[i][j] = *temp_coarse_flow;
							//swap the elements
					}
					k++;
						// look at next x-axis element in row
				}
				k = 1;
					// reset k because we successfuly ignored initial temp element
			}
		}
	}
	free(temp_coarse_flow);
	
	//#####VALGRIND ERROR CHECK#####
	//printf("currently on line 488\n");
	
	/* debug to check sort actually worked
	printf("scores matrix\n");
    for (int i = 1; i < (resolution+1); i++) {
		for (int j = 1; j < (resolution+1); j++) {
			printf("%f ",coarse_flows[i][j].score);
		}
		printf("\n");
	}*/
	
    // print all data to task2.csv
    // modified from ENGR3K3>L2>filo_io.c
    FILE* fp = safe_fopen("task2.csv", "w");
    if (fprintf(fp, "rho,u,v,x,y,S\n")<0){
    	perror("file write error");
    	exit(EXIT_FAILURE); 
	}   	
	
    for (int j = 1; j < (resolution+1); j++) {
		for (int i = 1; i < (resolution+1); i++) {
			if (fprintf(fp, "%f,%f,%f,%f,%f,%f\n", coarse_flows[i][j].rho, coarse_flows[i][j].u, coarse_flows[i][j].v, coarse_flows[i][j].x, coarse_flows[i][j].y, coarse_flows[i][j].score) < 0){
				perror("file write error");
				exit(EXIT_FAILURE);  
			}
		}    
    }	
        
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 517\n");
		// 1 valgrind error after this line 
	
	
    fclose(fp);
    
    /* clear the malloc for the 2D array */
	for (int i = 1; i < (resolution+1); i++) {
		/* some printf statements for debug */
		//printf("About to free coarse_flows[%d]\n",i);
        free(coarse_flows[i]);
        //printf("Free success\n");
        	//segmentation fault WAS happening here bc printf never runs. nothing was actually freed 
        	//fault has since been fixed
    }
	free(coarse_flows);
	
}

void searching(const char* flow_file)
{
    //printf("searching() - IMPLEMENT ME!\n");
    
    /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 1. malloc for 1D array of structs for domain centreline flow points
    //		1.1. read through file once, count how many y=0.0 there are
    //		1.2. close file
    //		1.3. open file and go to (2)
    // 2. loop that scanfs data from file
    // 		2.1. if point y=/=0.0, ignore it!
    //		2.2. directly store these points in 1D array of structs
    //		2.3. while doing the scanf, find out the highest rho*u
    //			##DONT NEED TO DO## Will find out max rho*u when doing the sort in (4) 
    // 3. close file
    // 4. sort 1D array in ascending order of rho*u
    //		4.1. note the highest rho*u
    //		4.2. set this as a variable and make sure it's ignored during searches later steps
    // 5. malloc linked list
    // 6. malloc binary search tree
    // 7. put the sorted data into a linked list
    // 8. put sorted 1D array somehow into BST perfectly balanced??
    // 		8.1. how to do this if the array is sorted, do we unsort or not allowed????// 8. 
    // 9. linear search array to find value closest to 40% * max rho*u
    //		9.1. write all checked values inc. found exact value to task3.csv
    //		9.2. output time taken for search//
    // 10. binary search array to find same value closest to 40% * max rho*u
    //		10.1. write all checked values inc. found exact value to task3.csv
    //		10.2. output time taken for search//
    // 11. linear search linked list to find same value closest to 40% * max rho*u
    //		11.1. write all checked values and inc. found exact value to task3.csv
    //		11.2. output time taken for search//
    // 12. search balanced BST to find same value closest to 40% * max rho*u
    //		12.1. write all checked values inc. found exact value to task3.csv
    //		12.2. output time taken for search//
    // 13. malloc free() 1D array, linked list, and BST
    // ###I slightly changed the order of some steps, but each individual step is still completed
    
  
    
    /* Open flow_data.csv file */
    FILE *f;
	f = safe_fopen(flow_file,"r");
    fscanf(f,"rho,u,v,x,y");
    	//ignore first line
    int y_zero_ctr = 0; //count amount of y=0.0 flow pts
    double temp_y=1.0;
    
    
    // read through file once, count how many y=0.0 there are
    for (int i=0;i<(TOTAL_FLOW_PTS);i++){
		fscanf(f,"%*f,%*f,%*f,%*f,%lf",&temp_y);
			//ignore all values except y
    	if (temp_y == 0.0){
    		y_zero_ctr += 1; 
    		//printf("i = %d, y = %lf\n",i,temp_y);
    	}	
    }
    //printf("num of flow points where y=0.0 --> %d\n",y_zero_ctr);
    	//the 209287th point in file (i=209286) is first y=0.0 flow point
    	//the 210343rd point in file (i=210342) is last y=0.0 flow pt
    	//in total 1057 points
    
    fclose(f);
    
   
    // malloc for 1D array of structs for domain centreline flow points
    flow_t3* flows = NULL;
	flows = (flow_t3*) malloc(y_zero_ctr*sizeof(flow_t3));
		// it's y_zero_ctr NOT y_zero_ctr+1 because I'm gonna try analysis from a[0] instead of a[1]
	assert(flows != NULL);
	

	//### Re-open file
	f = safe_fopen(flow_file,"r");
	if (fseek(f, 0, SEEK_SET)) {
		//reset file to beginning so I can read again
		perror("file reset error");
    	exit(EXIT_FAILURE); 
	}
	
	//printf("currently on line 643\n");
		
    fscanf(f,"rho,u,v,x,y");
    	//ignore first line
    	
    
    // loop that scanfs data from file
    double rho, u, v, x, y;
    int i_array = 0;
    	//to do storage of flows[i_array] from i_array = 0 to i_array = y_zero_ctr-1 
    	
    for (int i=0; i<(TOTAL_FLOW_PTS) ;i++){
    	fscanf(f,"%lf,%lf,%lf,%lf,%lf",&rho, &u, &v, &x, &y);
    	
    	//printf("currently on line 653, i = %d\n", i);
    	//printf("%f,%f,%f,%f,%f\n", rho, u, v, x, y);
    	if (y == 0.0){
    		// if point y=/=0.0, ignore it!

    		//printf("currently on line 667\n");
    		//printf("currently i = %d\n", i);
    		flows[i_array].rho = rho;
    		flows[i_array].u = u;
    		flows[i_array].v = v;
    		flows[i_array].x = x;
    		flows[i_array].y = y;
    		flows[i_array].rho_u = rho*u;
    		//printf("i_array = %d, rho*u = %0.8lf\n", i_array, flows[i_array].rho_u);
    		//printf("%f,%f,%f,%f,%f\n", flows[i_array].rho, flows[i_array].u, flows[i_array].v, flows[i_array].x, flows[i_array].y);
    		
    		//printf("element i_array = %d is processed\n", i_array);
    			//if all is correct, the final i_array value printed will be y_zero_ctr-1 aka 1057-1 = 1056
    			//all is correct!
    		i_array++;
    	}    
    }
    
    // close file
    fclose(f);
    
    
    //### sort 1D array in ascending order of rho*u
    qsort(flows, y_zero_ctr, sizeof(flow_t3), cmp_t3);
    	//qsort is like magic. I honestly don't know how one line can make something so simple, yet here I am
    
    // testing that it's actually properly sorted --> IT IS!!
    /*for(int i=0; i<y_zero_ctr; i++){
    	printf("i = %d, rho*u = %0.8lf\n", i, flows[i].rho_u);
    }*/
    
        //note the highest rho*u
    double rho_u_max = flows[y_zero_ctr-1].rho_u;
    	//max rho*u value
    	//make sure it's ignored during searches later steps
    //double forty_rho_u_max = rho_u_max*0.4;
    	//40% of max rho*u

    //printf("max rho*u = %0.8lf\n", rho_u_max);
    
    
    //### create a task3.csv. have it ready for writing
    FILE* fp = safe_fopen("task3.csv", "w"); 	   
	
    //### linear search array to find value closest to 40% * max rho*u
    //		9.2. output time taken for search
    int found_val=0, i=0;
    	//tracks whether we've found the val closest to 40% * rho_u_max
    struct timeval start;
	struct timeval stop;
		//for use in printing search times
	gettimeofday(&start, NULL);	
	
    double difference = fabs(flows[i].rho_u - 0.4*rho_u_max);
    	//seeing absolute value of difference between initial flows[i].rho_u and 40%*rho_u_max
    	//initialised as the first element in ascending order array
    double rho_u_max_40 = 0.4*rho_u_max;
    	// 40% of rho_u_max
    double found_rho_u;
    	
    fprintf(fp,"%f", flows[i].rho_u);
    	//print out first element
    	
    //printf("currently on line 1034, rho_u_max = %f, ", rho_u_max);
	//printf("0.4*rho_u_max = %f\n", rho_u_max_40);
    i=1;
    	//start searching flows[i] from i=1 (second element)
    	//LINEAR SEARCH. goes through every element until we find desired one
    while (found_val==0){
    	if (fabs(flows[i].rho_u - rho_u_max_40)<=difference){
    		//if smaller than difference of previous element, it is closer to desired val
    		//just in case difference is equal, carry on search
    		difference = fabs(flows[i].rho_u - rho_u_max_40);
    		fprintf(fp,",%f", flows[i].rho_u);
    			
    	} else{
    		//means difference of curent element is greater than difference of previous element
    		//when difference begins to increase, we are moving past the rho*u we want
    		found_val=1;
    			//end loop
    		//printf("element before found element --> rho*u = %f\n", flows[i-2].rho_u);
    		//printf("found element --> rho*u = %f\n", flows[i-1].rho_u);
    		//printf("element after found element --> rho*u = %f\n", flows[i].rho_u);
    			//found element is indeed the closest = rho*u = 0.415171
    		found_rho_u = flows[i-1].rho_u;
    	}
    	i++;
    }
    fprintf(fp,"\n");
	gettimeofday(&stop, NULL);
		//linear search on array completed
	double elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0 * 1000.0;
		//note this is in microseconds, not milli
	elapsed_ms += (stop.tv_usec - start.tv_usec);
	printf("TASK 3 Array Linear Search:  %.2f microseconds\n", elapsed_ms);
    
	
	//#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1030\n");
    
       	
    /* modified binary search from ENGR3K3 github 'binary_search.c' */
    //#### binary search array to find same value closest to 40% * max rho*u
    int lo = 0;
    int hi = y_zero_ctr-1;
    gettimeofday(&start, NULL);
    
    while (lo < hi) {
        int m = (lo + hi) / 2;
        	//m is for middle basically. splitting the array in half sort of
        if (found_rho_u < flows[m].rho_u) {
            fprintf(fp,"%f,", flows[m].rho_u);
            hi = m;
            	//cut off right side of sorted array
        } else if (found_rho_u > flows[m].rho_u) {
            lo = m + 1;
            	//cut off left side of sorted array
            fprintf(fp,"%f,", flows[m].rho_u);
        } else {
        	//found found_rho_u aka closest rho_u to 40%*rho_u_max
            fprintf(fp,"%f\n", flows[m].rho_u);
            lo=hi;
            	//end loop
        }
    }
    
    gettimeofday(&stop, NULL);
		//binary search on array completed
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0 * 1000.0;
		//note this is in microseconds, not milli
	elapsed_ms += (stop.tv_usec - start.tv_usec);
	printf("TASK 3 Array Binary Search:  %.2f microseconds\n", elapsed_ms);
    
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1067\n");
	
	
    //#### create malloc linked list
    /* from workshop 2 */
    /* create new list object with regular 'free' as the del function */
    list_t* list = list_new(free);
    
    // put the sorted data into a linked list
    // modified from workshop 2
    /* insert some elements at the front */
    assert(list->num_elements == 0);
    for (int i = y_zero_ctr-1; i >= 0; i--) {
    	//list_push_front() continously puts next node at the beginning. means linked list becomes descending order
    	//to do the opposite, 'int i = 0; i < y_zero_ctr; i++' changed to 'int i = y_zero_ctr-1; i >= 0; i--'
        flow_t3* new_dat = (flow_t3*)malloc(sizeof(flow_t3));
        *new_dat = flows[i];
        	// to assign to a pointer of single struct, either do '*new_dat = flows[i]' OR
        	// 'new_dat->rho = flows[i].rho'
        	// 'new_dat->u = flows[i].u'
        	// etc etc
        	// ###ONLY when malloc is done to make array of structs can you use __[i].__ otherwise it's __->__
        	
        //printf("linked list node i=%d, rho*u = %0.8lf\n", i, new_dat->rho_u);
        list_push_front(list, new_dat);
    }
    assert(list->num_elements == y_zero_ctr);
    	//check all nodes are added successfully
    	
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1097\n");
    	
    //### linear search linked list to find same value closest to 40% * max rho*u
    gettimeofday(&start, NULL);
    
    /* conduct linear search by processing elements in the linked list */
    /* modified from workshop 2 */
    node_t* p = list->head;
    	//declare pointer to a node in list. set it to whatever node the head points to
    flow_t3* data_tmp = (flow_t3*) p->data;
    	//data is cast as a flow_t3 ptr because in node_t, it was initially defined as a void* ptr
    while (p != NULL) {
        if (data_tmp->rho_u == found_rho_u){
        	fprintf(fp,"%f\n", data_tmp->rho_u);
        		//FOUND!!
        		//write to file
        	p = NULL;
        } else{
        	fprintf(fp,"%f,", data_tmp->rho_u);
        		//write to file
        	p = p->next;
        	data_tmp = (flow_t3*) p->data;
        }
        	//move to next node in linked list
    }
    
    /* placeholder in case I need it. from workshop 2
    node_t* cur = list->head;
    for (size_t i = 0; i < 5; i++) {
        cur = cur->next;
    }*/
        
    gettimeofday(&stop, NULL);
		//linear search on linked list completed
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0 * 1000.0;
		//note this is in microseconds, not milli
	elapsed_ms += (stop.tv_usec - start.tv_usec);
	printf("TASK 3 List Linear Search:  %.2f microseconds\n", elapsed_ms);
    //fprintf(fp,"\n");
    
    //### free linked list
    list_free(list);  
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1141\n");
    
    //#### create a binary search tree bst
    bst_t* bst = bst_new(no_free, cmp_t3); // memory is held by the array itself
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1148\n");
    
    /* modified from workshop 3 */
    /* perfect insert */
    perfect_insert(bst, flows, 0, y_zero_ctr - 1);
    	// put sorted 1D array into BST. perfectly balanced
    //printf("num_elements = %d\n", bst->num_elements);
    assert(bst->num_elements == y_zero_ctr);
    
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1160\n");
    
    gettimeofday(&start, NULL);
    
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1165\n");    
    
	node_bst_t* tmp = bst->root;
		//declare pointer to a node_bst_t in BST. set it to whatever node the root points to
	data_tmp = (flow_t3*) tmp->data;
		//already defined data_tmp before. data type remains the same
    	//data is (or previously was) cast as a flow_t3 ptr because in node_bst_t (and node_t for list), it was initially defined as a void* ptr
 
    //#####VALGRIND ERROR CHECK#####
	//printf("currently on line 1174\n");
    	
    //### search balanced BST to find same value closest to 40% * max rho*u
    // modified from ENGR3K3 github bst.c. Also in L3	
	while(tmp != NULL) {
		if(found_rho_u < data_tmp->rho_u){
			//if our desired search val is less than current leaf/node in tree, move left
			fprintf(fp,"%f,", data_tmp->rho_u);
			tmp = tmp->left;
			data_tmp = (flow_t3*) tmp->data;
		} else if(data_tmp->rho_u < found_rho_u) {
			//if our desired search val is larger than current leaf/node in tree, move right
			fprintf(fp,"%f,", data_tmp->rho_u);
			tmp = tmp->right;
			data_tmp = (flow_t3*) tmp->data;
		} else {
			/* FOUND! */
			fprintf(fp,"%f\n", data_tmp->rho_u);
			tmp = NULL;
		}
		//move onto next node, unless desired node found, in which case while() ends
	}
    
    gettimeofday(&stop, NULL);
		//linear search on linked list completed
	elapsed_ms = (stop.tv_sec - start.tv_sec) * 1000.0 * 1000.0;
		//note this is in microseconds, not milli
	elapsed_ms += (stop.tv_usec - start.tv_usec);
	printf("TASK 3 BST Search:  %.2f microseconds\n", elapsed_ms);
    
    //### malloc free() BST
    bst_free(bst);
    
    // free 1d array. not needed anymore
    free(flows);
    	//the BST construction is made such that it's linked to flows even after the BSTs own creation
    	//if free(flows) is done after the creation of BST but before the search, there are valgrind errors despite program running perfectly
    	//solution is to free(flows) basically at the end

    fclose(fp);
    	//close task3.csv
}

void vortcalc(const char* flow_file)
{
    //printf("vortcalc() - IMPLEMENT ME!\n");
    
        /* MY LOGIC FOR SOLUTION */
    /* ~~~~~~~~~~~~~~~~~~~~~ */
    // 0. define a new flows_t4 struct. has an omega 'w' double element = vorticity
    //		0.1. ### dont think this is necessary. As loop (6) runs, just assign each flow pt to whatever (6.1) case
    //		0.2. (0.1) means I can use flow_t1 struct.  
    //		0.3. I might just add an omega 'w' vorticity element in a new flow_t4 struct just in case I get marked down otherwise
    // 1. open file
    //		1.1. read through file once, find out the following numbers
    //			1.1.2. find out number of unique y-values (maybe = 419629/1057 = 397?)
    //			1.1.3. find out how many flow points for each unique y-val (=1057)
    //				1.1.3.1. note the same x-values exist for a given y.
    //		1.2. close file
    // 2. malloc for 2D array of structs for flow points flows[i][j]. 2d array of structs is allowed as per A1 pdf
    //		2.1. for a given j (y), i increases--> x increases
    //		2.2. for a given i (x), j increases--> y increases
    //		2.3. result is an actual grid.
    // 3. open file again
    // 4. loop that fscanfs data from file into malloc'd 2D array
    // 		4.1. take initial y-val (=-20) for j=0, run through loop for scanf 1057 times for i=0 to 1057-1
    //		4.2. increment j
    //		4.3. repeat (4.1) with the next y-val and continue until j has run from j=0 to j=(1.1.2)=397-1
    // 5. close file
    // 6. loop through entire array and calculate vorticities noting the special cases in ass1 pdf
    // 		6.1. have a check and counters for abs|w| vorticities that lie in 0<=|w|<5, 5<=|w|<10, 10<=|w|<15, 15<=|w|<20, 20<=|w|<25 
    // 7. fopen (create with write mode w) task4.csv file
    //		7.1. print the counters (6.1) to file in format specified in ass1 pdf
    // 8. malloc free()
    // 9. bob's your uncle
    
    
    //## open file
    FILE *f = safe_fopen(flow_file,"r");
   
    fscanf(f,"rho,u,v,x,y");
    	//ignore first line
    	
    //##read through file, find out the following numbers
    int unique_y_ctr = 0; //count amount 
    	//number of unique y-values (maybe = 419629/1057 = 397?)
    	//AKA how many flow points exist for given x-val. AKA Array[i][j], j=0 to unique_y_ctr-1. AKA 'm' if grid is n x m
    int y_repeats_ctr = 0;
    	//how many flow points exist for each unique y-val (=1057). 
    	//AKA number of unique x-values. AKA Array[i][j], i=0 to y_repeats_ctr-1. AKA 'n' if grid is n x m 
    	//note the same x-values exist for a given y. it's a proper grid
    double temp_y_prev;
    fscanf(f,"%*f,%*f,%*f,%*f,%lf",&temp_y_prev);
    	//initialise second line to begin the comparisons
    double temp_y_current;
    unique_y_ctr++;
    y_repeats_ctr++;
    
    // read through file
    for (int i=1;i<(TOTAL_FLOW_PTS);i++){
    	//starting from i=1 because we've already scanned one flow point
		fscanf(f,"%*f,%*f,%*f,%*f,%lf",&temp_y_current);
			//ignore all values except y
    	if (temp_y_current == temp_y_prev){
    		//if another flow point for same given y
    		y_repeats_ctr++;
    	} else {
    		//if they're =/=, then we've moved onto next sequence of flow points for a given y
    		y_repeats_ctr=0;
    		y_repeats_ctr++;
    			//reset repeats counter then increment to include temp_y_current flow point
    		unique_y_ctr++;
    	}
    	temp_y_prev = temp_y_current;
    }
    //printf("there are %d unique y values and %d flow points for a given y\n",unique_y_ctr,y_repeats_ctr);
    
    //close file
    fclose(f);
    
    //### malloc for 2D array of structs for flow points flows[i][j]
	flow_t4** flows = NULL;
	flows = (flow_t4**) malloc((y_repeats_ctr)*sizeof(flow_t4));
	assert(flows != NULL);	
	for (int i = 0; i < (y_repeats_ctr); i++){
		flows[i] = (flow_t4*) calloc((unique_y_ctr),sizeof(flow_t4));
			// ####### it's resolution+1 because c starts arrays from a[0] but I want to start analysis from a[1]
			// since i've doing a bunch of += and /= later in some loops, i use calloc to initialise to 0
			// note that instead of malloc((resolution+1)*sizeof(coarse_flow_t2)), calloc is written as above w/ the comma
	}
    
    //### Open file again, reset position
    f = safe_fopen(flow_file,"r");
		//reopen file
	if (fseek(f, 0, SEEK_SET)) {
		//reset file to beginning so I can read again
		perror("file reset error");
    	exit(EXIT_FAILURE); 
	}
	fscanf(f,"rho,u,v,x,y");
    	//ignore first line
    
    //###loop that fscanfs data from file into malloc'd 2D array
    for(int j=0; j < unique_y_ctr; j++){
    	//run through all given y values, of which there are unique_y_ctr amount (=397)
    	for(int i=0; i < y_repeats_ctr; i++){
    		//for a given j (y), go through every ascending i (x)
    		//There are y_repeats_ctr (=1057) i elements for every j
    		fscanf(f,"%lf,%lf,%lf,%lf,%lf",&flows[i][j].rho, &flows[i][j].u, &flows[i][j].v, &flows[i][j].x, &flows[i][j].y);
    		
    		//printf("for flows[i][j]=flows[%d][%d], x = %f, y = %f\n", i, j, flows[i][j].x, flows[i][j].y);
    		//all data read successfully!!
    	}
    }
    
    //###close file
    fclose(f);
    
    int ctr_0_5=0;
    int ctr_5_10=0;
    int ctr_10_15=0;
    int ctr_15_20=0;
    int ctr_20_25=0;
    	//counters for abs|w| vorticities that lie in 0<=|w|<5, 5<=|w|<10, 10<=|w|<15, 15<=|w|<20, 20<=|w|<25 
    
    	
    //### loop through entire array and calculate vorticities noting the special cases
    for(int j=0; j < (unique_y_ctr); j++){
    	//run through all given y values for j = 0:m-1 = 0:(unique_y_ctr-1) = 0:(397-1)
    	for(int i=0; i < (y_repeats_ctr); i++){
    		//for a given j (y), go through every i (x) = 0:n-1 = 0:(y_repeats_ctr-1) = 0:(1057-1)
    		
    		if(i==(y_repeats_ctr-1) && j!=(unique_y_ctr-1)){
    			//1st special case
    			flows[i][j].w = (flows[i][j].v - flows[i-1][j].v)/(flows[i][j].x - flows[i-1][j].x) - (flows[i][j+1].u - flows[i][j].u)/(flows[i][j+1].y - flows[i][j].y);
    		
    		}else if (i!=(y_repeats_ctr-1) && j==(unique_y_ctr-1)){
    			//2nd special case
    			flows[i][j].w = (flows[i+1][j].v - flows[i][j].v)/(flows[i+1][j].x - flows[i][j].x) - (flows[i][j].u - flows[i][j-1].u)/(flows[i][j].y - flows[i][j-1].y);
    		
    		}else if (i==(y_repeats_ctr-1) && j==(unique_y_ctr-1)){
    			//3rd special case
    			flows[i][j].w = (flows[i][j].v - flows[i-1][j].v)/(flows[i][j].x - flows[i-1][j].x) - (flows[i][j].u - flows[i][j-1].u)/(flows[i][j].y - flows[i][j-1].y);
    		
    		}else{
    			//not special case
				flows[i][j].w = (flows[i+1][j].v - flows[i][j].v)/(flows[i+1][j].x - flows[i][j].x) - (flows[i][j+1].u - flows[i][j].u)/(flows[i][j+1].y - flows[i][j].y);
					//calculate vorticity using general formula from A1 pdf
    		}
    		
			//check if absolute value of vorticity lies within the ranges we're interested in
			if (0 <= fabs(flows[i][j].w) && fabs(flows[i][j].w) < 5){
				ctr_0_5++;
			}else if(5 <= fabs(flows[i][j].w) && fabs(flows[i][j].w) < 10){
				ctr_5_10++;
			}else if(10 <= fabs(flows[i][j].w) && fabs(flows[i][j].w) < 15){
				ctr_10_15++;
			}else if(15 <= fabs(flows[i][j].w) && fabs(flows[i][j].w) < 20){
				ctr_15_20++;
			}else if(20 <= fabs(flows[i][j].w) && fabs(flows[i][j].w) < 25){
				ctr_20_25++;
			}
    	}
    }
    
    
    //###fopen (create with write mode w) task4.csv file
    //		7.1. print the counters (6.1) to file in format specified in ass1 pdf
    FILE* fp = safe_fopen("task4.csv", "w");
    if (fprintf(fp, "threshold,points\n")<0){
    	perror("file write error");
    	exit(EXIT_FAILURE); 
	}  
	
	//##############DO PRINTING###############
	fprintf(fp, "5,%d\n", ctr_0_5);
	fprintf(fp, "10,%d\n", ctr_5_10);
	fprintf(fp, "15,%d\n", ctr_10_15);
	fprintf(fp, "20,%d\n", ctr_15_20);
	fprintf(fp, "25,%d\n", ctr_20_25);

    
	//### close file
	fclose(fp);
	
    // clear the malloc for the 2D array
	for (int i = 0; i < (y_repeats_ctr); i++) {
		//printf("About to free coarse_flows[%d]\n",i);
        free(flows[i]);
        	//each flows[i] corresponds to a bunch of flows[i][j], where j=0 to unique_y_ctr-1
        	//i goes from i = 0 to i = y_repeats_ctr-1
        //printf("Free success\n");
    }
	free(flows);

}
/* and now my watch has ended. thanks for reading!! */

/* ~~~~~~~~~~~ fin ~~~~~~~~~~ */
