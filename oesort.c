// Sample MPI Program for Parallel Programming 2015 Lab1
// You can reuse this code as a template for homework 1.
// Or you can write your own from scratch!
// It's up to you.

#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>

#define ROOT 0
#define is_qs 100
inline void printall(int *array, int length){	
	int i;
	for(i=0;i<length;i++){
		printf("%d ",array[i]);
	}
	putchar('\n');
}
int cmp(const void* a, const void* b){
	if(*(int*)a > *(int*)b)return 1;
	else if(*(int*)a < *(int*)b)return -1;
	return 0; 
}
void insertionsort(int* array, int length)
{
	int n, j, i;
	for (i = 1; i < length; ++i){
		n = array[i];
		for (j = i - 1; j >= 0 && array[j] > n; --j)
			array[j + 1] = array[j];
		array[j + 1] = n;
	}
}

void qsort_int(int* array, int length){
	qsort((void*)array, length, sizeof(int), cmp);
}
int main (int argc, char *argv[]) {
	int rank, size;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Status status;
	if (argc < 4) {
		if (rank == ROOT) {
			fprintf(stderr, "Insufficient args\n");
			fprintf(stderr, "Usage: %s N input_file", argv[0]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return 0;
	}

	int N = atoi(argv[1]), alloc_num;
	const char *inName = argv[2];
	const char *outName = argv[3];
	double start, finish;
	int *root_arr; // for root process (which rank == 0) only

	// Part 1: Read file
	/* Note: You should deal with cases where (N < size) in Homework 1 */
	int rc;
	MPI_File fp;
	rc = MPI_File_open(MPI_COMM_WORLD, inName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp); 
	if(rc != MPI_SUCCESS){
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Offset total_number_of_bytes;
	MPI_File_get_size(fp, &total_number_of_bytes);
	if(total_number_of_bytes/sizeof(int)<N){
		puts("N is bigger than testcase in input file, read to the end");
		N = total_number_of_bytes/sizeof(int);
	}
	// sheu if N < # of processes?
	int *array;
	if(N < size){
		if(rank != ROOT){
			MPI_Finalize();
			exit(0);
		}
		else{
			// sheu todo
			puts("too less input, use insertion sort by ROOT process");
			array = (int*)malloc(sizeof(int)*N);
			MPI_Status status;
			MPI_File_read(fp, array, N, MPI_INT, &status);
			insertionsort(array,N);
			printall(array,N);
			MPI_Finalize();
		}
	}
	else{
		alloc_num = N / size;
		if(N%size==0){
			//todo
			MPI_File_seek(fp, (MPI_Offset)rank*alloc_num, MPI_SEEK_SET);
			array = (int*)malloc(sizeof(int)*alloc_num);
			start = MPI_Wtime();
			MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			printf("rank: %2d io time: %lf\n",rank,finish-start);
		}
		else{
			alloc_num+=1;
			if(rank*alloc_num>=N){
				printf("rank: %2d others take my job, exit\n",rank);
				MPI_Finalize();
				exit(0);
			}
			MPI_File_seek(fp, (MPI_Offset)rank*alloc_num*sizeof(int), MPI_SEEK_SET);
			array = (int*)malloc(sizeof(int)*alloc_num);
			if(rank*alloc_num+alloc_num>N){
				alloc_num = N-rank*alloc_num;
				start = MPI_Wtime();
				MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
				finish = MPI_Wtime();
				printf("rank: %2d io time: %lf\n",rank,finish-start);
			}
			else{
				start = MPI_Wtime();
				MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
				finish = MPI_Wtime();
				printf("rank: %2d io time: %lf\n",rank,finish-start);
			}
			//test
			insertionsort(array,alloc_num);
			printall(array,alloc_num);
			MPI_Finalize();
			exit(0);
		}
		
		
	}
	MPI_Finalize();
	exit(0);

/*	if (rank == ROOT) {
		struct stat st;
		int fd = 0;
		fd = open(argv[2],O_RDONLY);
		fstat(fd, &st);
		if(fd != -1){
			puts("N is bigger than testcase in input file, read to the end");
		N = st.sts_size/sizeof(int);
		}
		else{
			puts("ferror QQ");
			MPI_Barrier(MPI_COMM_WORLD);
			MPI_Finalize();
		}
		close(fd);
	}
	MPI_Barrier(MPI_COMM_WORLD);
*/
/*		root_arr = malloc(N * sizeof(int));
		FILE *f = fopen(inName, "rb"); // sheu may need change
		fread(root_arr, sizeof(int), N, f);
		fclose(f);
	}
*/

	// Part 2: Scatter
	/* Note: You should deal with cases where (N % size != 0) in Homework 1 */
	int num_per_node = N / size;
	int *local_arr = malloc (num_per_node * sizeof(int));
	// Ref:
	// MPI_Scatter (send_buf, send_count, send_type, recv_buf, recv_count, recv_type, root, comm)
	MPI_Scatter    (root_arr, num_per_node, MPI_INT, local_arr, num_per_node, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (rank == ROOT) {
		free (root_arr);
	}

	// Part 3: Calculate the sum of local_arr
	int sum = 0;
	int i;
	for (i = 0; i < num_per_node; i++) {
		sum += local_arr[i];
	}

	free (local_arr);

	// Part 4: Accumulate the result to root process
	int accumulated;
	// Ref:
	// MPI_Reduce (send_buf, recv_buf, count, data_type, op, root, comm)
	MPI_Reduce    (&sum, &accumulated, 1, MPI_INT, MPI_SUM, ROOT, MPI_COMM_WORLD);

	if (rank == ROOT) {
		printf("SUM: %d\n", accumulated);
	}

	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
