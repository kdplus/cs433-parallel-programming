#include <stdio.h>
#include <string.h>
#include <mpi.h>

const int MAX_STRING = 100;
#define comm MPI_COMM_WORLD

const int INFINITY = 1000000;

void Read_matrix(int mat[], int n);
void Print_matrix(int mat[], int n);
void Print_dists(int dist[], int n);
void Print_paths(int pred[], int n);
int  Find_min_dist(int dist[], int known[], int n);
void Dijkstra_Para(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int my_rank, int N);

int main(void) {
	char greeting[MAX_STRING];
	int comm_sz;
	int my_rank;
	
	int  n, loc_n;
	int *mat, *dist, *pred;
	int *loc_mat, *loc_dist, *loc_pred;


	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);


	if (my_rank != 0) {
	/*
		sprintf(greeting, "Greetings from process %d of %d!", my_rank, comm_sz);
		MPI_Send(greeting, strlen(greeting) + 1, MPI_CHAR, 0, 0, MPI_COMM_WORLD);
	*/
		MPI_Recv(&n, 1, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	// Receive the size of matrix from process 0
		loc_n = n/comm_sz;	// Divide computation into comm_sz parts
		loc_mat = malloc(n*n/loc_n*sizeof(int));
		loc_dist = malloc(n/loc_n*sizeof(int));
		loc_pred = malloc(n/loc_n*sizeof(int));
		MPI_Scatter(mat, loc_n*n, MPI_INT, loc_mat, loc_n*n, MPI_INT, 0, MPI_COMM_WORLD);	// Scatter matrix mat whose size is loc_n * n
		Dijkstra_Para(loc_mat, loc_dist, loc_pred, loc_n, my_rank, n); 
		MPI_Gather(loc_dist, loc_n, MPI_INT, dist, loc_n, MPI_INT, 0, MPI_COMM_WORLD); // Form the complete dist
		MPI_Gather(loc_pred, loc_n, MPI_INT, pred, loc_n, MPI_INT, 0, MPI_COMM_WORLD); // Form the complete pred
	} else {
	/*	
		printf("Greetings from process %d of %d!\n", my_rank, comm_sz);
		for(int q = 1; q < comm_sz; q++){
			MPI_Recv(greeting, MAX_STRING, MPI_CHAR, q, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			printf("%s\n", greeting);
		}
	*/	
		printf("How many vertices?\n");
		scanf("%d", &n);
		for(int i = 1; i < comm_sz; ++i )
			MPI_Send(&n, 1, MPI_INT, i, 0, MPI_COMM_WORLD);
		mat = malloc(n*n*sizeof(int));
		dist = malloc(n*sizeof(int));
		pred = malloc(n*sizeof(int));
		loc_n = n/comm_sz;
		loc_mat = malloc(n*n/loc_n*sizeof(int));
		loc_dist = malloc(n/loc_n*sizeof(int));
		loc_pred = malloc(n/loc_n*sizeof(int));
		
		printf("Enter the matrix\n");
		Read_matrix(mat, n);

		MPI_Scatter(mat, loc_n*n, MPI_INT, loc_mat, loc_n*n, MPI_INT, 0,
				MPI_COMM_WORLD); // Scatter matrix mat whose size is loc_n * n
		Dijkstra_Para(loc_mat, loc_dist, loc_pred, loc_n, my_rank, n);

		MPI_Gather(loc_dist, loc_n, MPI_INT, dist, loc_n, MPI_INT, 0,
				MPI_COMM_WORLD);  // Form the complete dist
		MPI_Gather(loc_pred, loc_n, MPI_INT, pred, loc_n, MPI_INT, 0,
				MPI_COMM_WORLD);  // Form the complete pred
		printf("The distance from 0 to each vertex is:\n");
		Print_dists(dist, n);
		printf("The shortest path from 0 to each vertex is:\n");
		Print_paths(pred, n);

		free(mat);
		free(dist);
		free(pred);
	}
	MPI_Finalize();
	return 0;
}


/*-------------------------------------------------------------------
 * Function:  Read_matrix
 * Purpose:   Read in the adjacency matrix
 * In arg:    n
 * Out arg:   mat
 */
void Read_matrix(int mat[], int n) {
   int i, j;

   for (i = 0; i < n; i++){
      for (j = 0; j < n; j++){
		  scanf("%d", &mat[i*n+j]);
	  }
   }
}  /* Read_matrix */

/*-------------------------------------------------------------------
 * Function:  Print_matrix
 * Purpose:   Print the contents of the matrix
 * In args:   mat, n
 */
void Print_matrix(int mat[], int n) {
   int i, j;

   for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++)
         if (mat[i*n+j] == INFINITY)
            printf("i ");
         else
            printf("%d ", mat[i*n+j]);
      printf("\n");
   }
}  /* Print_matrix */

/*-------------------------------------------------------------------
 * Function:    Dijkstra_Para
 * Purpose:     Apply Dijkstra_Para's algorithm to the matrix loc_mat
 * In args:     loc_n:  the number of vertices
 *              loc_mat:  adjacency matrix for the graph
 * Out args:    loc_dist:  loc_dist[v] = distance 0 to v.
 *              loc_pred:  loc_pred[v] = predecessor of v on a 
 *                  shortest path 0->v.
 */
void Dijkstra_Para(int loc_mat[], int loc_dist[], int loc_pred[], int loc_n, int my_rank, int N) {
	int i, loc_u, v, *known, new_dist;
	int my_min[2], glbl_min[2];
	/* known[v] = true, if the shortest path 0->v is known */
	/* known[v] = false, otherwise                         */
	known = malloc(loc_n*sizeof(int));
	/* Initialize d and p */
	loc_dist[0] = 0; loc_pred[0] = 0; known[0] = 1; 
	for (v = 0; v < loc_n; v++) { // In paralleliztion, v should start from 0, not 1
		loc_dist[v] = loc_mat[0 + v*N];
		loc_pred[v] = 0;
		known[v] = 0;
	}

#     ifdef DEBUG
	  printf("i = 0\n");
	  Print_dists(loc_dist, loc_n);
#     endif

	/* On each pass find an additional vertex */
	/* whose distance to 0 is known           */
	for (i = 1; i < loc_n; i++) {
	  	loc_u = Find_min_dist(loc_dist, known, loc_n);
		int global_u = loc_u + my_rank * loc_n;	// Find global vertex number
		my_min[0] = loc_dist[loc_u]; // Distance from 0 to loc_u
		my_min[1] = global_u; // Store global vertex number
		
		MPI_Allreduce(my_min, glbl_min, 1, MPI_2INT, MPI_MINLOC,
				MPI_COMM_WORLD); // glbl_min[0]: current global minimum dist, glbl_min[1]: global vertex number of this chosen vertex
	  	if (global_u == glbl_min[1]) known[loc_u] = 1; // Analyze whether the chosen vertex is in this process

		for (v = 0; v < loc_n; v++) 
			if (!known[v]) {
				new_dist = glbl_min[0] + loc_mat[glbl_min[1] + v*N]; // u + v * N
			/*
				printf("________________________\n");
				printf("%d\n", glbl_min[1]);
				printf("%d----%d\n", loc_dist[glbl_min[1]], loc_mat[glbl_min[1]*n+v]);
				printf("________________________\n");
			*/
				if (new_dist < loc_dist[v]) {
					loc_dist[v] = new_dist;
					loc_pred[v] = glbl_min[1];
				}
			}

#     ifdef DEBUG
	  printf("i = %d\n", i);
	  Print_dists(loc_dist, loc_n);
#     endif
	} /* for i */

	free(known);
}  /* Dijkstra_Para */

/*-------------------------------------------------------------------
 * Function:    Find_min_dist
 * Purpose:     Find the vertex u with minimum distance to 0
 *              (dist[u]) among the vertices whose distance 
 *              to 0 is not known.
 * In args:     dist:  dist[v] = current estimate of distance
 *                 0->v
 *              known:  whether the minimum distance 0-> is
 *                 known
 *              n:  the total number of vertices
 * Ret val:     The vertex u whose distance to 0, dist[u]
 *              is a minimum among vertices whose distance
 *              to 0 is not known.
 */
int Find_min_dist(int dist[], int known[], int n) {
   int v, u, best_so_far = INFINITY;

   for (v = 1; v < n; v++)
      if (!known[v])
         if (dist[v] < best_so_far) {
            u = v;
            best_so_far = dist[v];
         }

   return u;
}  /* Find_min_dist */


/*-------------------------------------------------------------------
 * Function:    Print_dists
 * Purpose:     Print the length of the shortest path from 0 to each
 *              vertex
 * In args:     n:  the number of vertices
 *              dist:  distances from 0 to each vertex v:  dist[v]
 *                 is the length of the shortest path 0->v
 */
void Print_dists(int dist[], int n) {
   int v;

   printf("  v    dist 0->v\n");
   printf("----   ---------\n");
                  
   for (v = 1; v < n; v++)
      printf("%3d       %4d\n", v, dist[v]);
   printf("\n");
} /* Print_dists */  


/*-------------------------------------------------------------------
 * Function:    Print_paths
 * Purpose:     Print the shortest path from 0 to each vertex
 * In args:     n:  the number of vertices
 *              pred:  list of predecessors:  pred[v] = u if
 *                 u precedes v on the shortest path 0->v
 */
void Print_paths(int pred[], int n) {
   int v, w, *path, count, i;

   path =  malloc(n*sizeof(int));

   printf("  v     Path 0->v\n");
   printf("----    ---------\n");
   for (v = 1; v < n; v++) {
      printf("%3d:    ", v);
      count = 0;
      w = v;
      while (w != 0) {
         path[count] = w;
         count++;
         w = pred[w];
      }
      printf("0 ");
      for (i = count-1; i >= 0; i--)
         printf("%d ", path[i]);
      printf("\n");
   }

   free(path);
}  /* Print_paths */
