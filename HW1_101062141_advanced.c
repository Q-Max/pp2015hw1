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
#include <string.h> // memcmp
#define ROOT 0
#define NOTSORTED 444
#define IS_QS 100
inline void printall(int *array, int length){	
	int i;
	for(i=0;i<length;i++){
		printf("%d ",array[i]);
	}
	putchar('\n');
}
inline void swap(int* a, int* b){
	*a = *a ^ *b;
	*b = *a ^ *b;
	*a = *a ^ *b;
}
int cmp(const void* a, const void* b){
	if(*(int*)a > *(int*)b)return 1;
	else if(*(int*)a < *(int*)b)return -1;
	return 0; 
}
inline void insertionsort(int* array, int length)
{
	int n, j, i;
	for (i = 1; i < length; ++i){
		n = array[i];
		for (j = i - 1; j >= 0 && array[j] > n; --j)
			array[j + 1] = array[j];
		array[j + 1] = n;
	}
}

inline void qsort_int(int* array, int length){
	qsort((void*)array, length, sizeof(int), cmp);
}
void singleOESort(int* array, int length){
	int i,sorted=0;
	while(!sorted){
		sorted = 1;
		for(i=0;i+1<length;i+=2){
			if(array[i]>array[i+1]){
				swap(&array[i],&array[i+1]);
				sorted = 0;
			}
		}
		for(i=1;i+1<length;i+=2){
			if(array[i]>array[i+1]){
				swap(&array[i],&array[i+1]);
				sorted = 0;
			}
		}
	}
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
			fprintf(stderr, "Usage: %s N input_file output_file", argv[0]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return 0;
	}

	int N = atoi(argv[1]), alloc_num, former_alloc_num=0, last_alloc_num=0;
	const char *inName = argv[2];
	const char *outName = argv[3];
	double start, finish, iotime = 0, commtime = 0, io_all, comm_all;
	int *root_ptr; // for root process (which rank == 0) only
  
	// Part 1: Read file
	/* Note: You should deal with cases where (N < size) in Homework 1 */
	int rc, i, trend = NOTSORTED;
	int *array;
	MPI_File fp;
	MPI_File fh;
	MPI_Offset my_offset;
	rc = MPI_File_open(MPI_COMM_WORLD, inName, MPI_MODE_RDONLY, MPI_INFO_NULL, &fp); 
	if(rc != MPI_SUCCESS){
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Offset total_number_of_bytes;
	MPI_File_get_size(fp, &total_number_of_bytes);
	if(total_number_of_bytes/sizeof(int)<N){
		if(rank==ROOT)
			puts("N is bigger than testcase in input file, read to the end");
		N = total_number_of_bytes/sizeof(int);
	}
	// sheu
    // todo
    // detect whether the file is sorted
    //--------------------------------------------------------------------------------
	if(rank!=ROOT){
		MPI_File_open(MPI_COMM_WORLD, outName, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);		
	}
	else{
		alloc_num = N;
		MPI_File_seek(fp,(MPI_Offset)0, MPI_SEEK_SET);
		array = (int*)malloc(sizeof(int)*alloc_num);
		start = MPI_Wtime();
		MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
		finish = MPI_Wtime();
		iotime += finish-start;
		MPI_File_open(MPI_COMM_WORLD, outName, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);		
		if(N>=2){
#ifdef DEBUG
			printall(array, alloc_num);
#endif
			if(array[0]>array[1])
				trend = -1;
			else if(array[0]==array[1])
				trend = 0;
			else
				trend = 1;
			for(i=2;i<alloc_num;i++){
				if(array[i-1]>array[i]&&(trend==-1||trend==0)){
					trend = -1;
					continue;
				}
				else if(array[i-1]==array[i]&&trend==0)
					continue;
				else if(array[i-1]<array[i]&&(trend==1||trend==0)){
					trend = 1;
					continue;
				}
				else{
					trend = NOTSORTED;
					break;
				}
			}
#ifdef DEBUG
			printf("%d\n", trend);
#endif
			if(trend==1||trend==0){
				printf("sorted file\n");
				my_offset = 0;
				printall(array, alloc_num);
				start = MPI_Wtime();
				MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
				finish = MPI_Wtime();
				iotime += finish - start;
				printf("iotime   : %8.5lf\ncommtime : %8.5lf\n",iotime,commtime);			
			}
			else if(trend==-1){
				printf("descending sorted file");
				root_ptr = (int*)malloc(sizeof(int)*alloc_num);
				for(i=0;i<alloc_num;i++){
					root_ptr[i] = array[alloc_num-i-i];
				}
				my_offset = 0;
				start = MPI_Wtime();
				MPI_File_write_at(fh, my_offset, root_ptr, alloc_num, MPI_INT, &status);
				finish = MPI_Wtime();
				iotime += finish - start;
				printf("iotime   : %8.5lf\ncommtime : %8.5lf\n",iotime,commtime);
			}
		}
		else{
			// N==1
			trend = 0;
			my_offset = 0;
			start = MPI_Wtime();
			MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			iotime += finish - start;
			printf("iotime/t:%8.5lf/ncommtime/t:%8.5lf\n",iotime,commtime);
			
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Bcast(&trend, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	if(trend==-1||trend==0||trend==1){
		MPI_Finalize();
		exit(0);
	}
	// read data in memory
	// sheu if N < 2 x  # of processes, root take over all computation
	
	alloc_num = N/size;
	if(2*size>N){
		// if N < 2 x size, root take over
		if(rank!=ROOT){			
			MPI_Finalize();
			exit(0);
		}
		else{
			alloc_num = N;
			MPI_File_seek(fp,(MPI_Offset)0, MPI_SEEK_SET);
			array = (int*)malloc(sizeof(int)*alloc_num);
			start = MPI_Wtime();
			MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			iotime = finish-start;
			singleOESort(array, alloc_num);
#ifdef DEBUG
			printall(array, alloc_num);      
#endif
			my_offset = 0;
			start = MPI_Wtime();
			MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			iotime += finish - start;
			printf("iotime/t:%8.5lf/ncommtime/t:%8.5lf\n",iotime,commtime);
			MPI_Finalize();
			exit(0);
		}
	}
/*	else if(alloc_num==1){
		// if N < 2*size, root take over
		if(rank!=ROOT){
			MPI_Finalize();
			exit(0);
		}
		else{
			alloc_num = N;
			MPI_File_seek(fp,(MPI_Offset)0, MPI_SEEK_SET);
			array = (int*)malloc(sizeof(int)*alloc_num);
			start = MPI_Wtime();
			MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			iotime = finish-start;
			singleOESort(array, alloc_num);
			//printall(array, alloc_num);
			my_offset = 0;
			start = MPI_Wtime();
			MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			iotime += finish - start;
			printf("iotime/t:%8.5lf/ncommtime/t:%8.5lf\n",iotime,commtime);
			MPI_Finalize();
			exit(0);
		}
	}*/
	else if((alloc_num)%2){
		// if alloc_num is odd number
		// every process alloc_num-- to guarantee every process has even elements
		// except last node
		alloc_num--;
		former_alloc_num = alloc_num;
		MPI_File_seek(fp,(MPI_Offset)alloc_num*rank*sizeof(int), MPI_SEEK_SET);
		if(rank==size-1){
			alloc_num = N - alloc_num * rank;
		}
		last_alloc_num = N - alloc_num * rank;
		array = (int*)malloc(sizeof(int)*alloc_num*2);
		start = MPI_Wtime();
		MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
		finish = MPI_Wtime();
		iotime += finish-start;
	}
	else{
		// last process need deal with more inputs
		former_alloc_num = alloc_num;
		MPI_File_seek(fp, (MPI_Offset)rank*alloc_num*sizeof(int), MPI_SEEK_SET);
		if(rank==size-1){			
			alloc_num = N - alloc_num * rank;
		}	
		last_alloc_num = N - alloc_num * rank;
		array = (int*)malloc(sizeof(int)*alloc_num*2);
		start = MPI_Wtime();
		MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
		finish = MPI_Wtime();
		iotime += finish-start;
	}
	
	// sheu
	// todo
	// sort and communicate with other
	//--------------------------------------------------------------------------------
	int tmp1,tmp2,sorted_temp;
	int sorted=0;
#ifdef DEBUG
	int *num_ptr, *pos_ptr; 
	if(rank==ROOT){
		root_ptr = (int*)calloc(0,sizeof(int)*N);	
		num_ptr = (int*)malloc(sizeof(int)*size);
		pos_ptr = (int*)malloc(sizeof(int)*size);
		for(i=0;i<size;i++){
			num_ptr[i] = former_alloc_num;
			pos_ptr[i] = i * former_alloc_num;
		}
		num_ptr[size-1] = last_alloc_num;
	}
#endif
	int *temp_array = malloc(sizeof(int)*alloc_num*2);
	int *sorted_array = malloc(sizeof(int)*alloc_num*2);
	if(alloc_num>IS_QS)
		quicksort=1;
	if(quicksort)
		qsort_int(array,alloc_num);
	else
		insertionsort(array,alloc_num);
	while(!sorted){
		sorted=1;    
		//printf("rank: %2d %d %d\n", rank, array[0], array[1]);	
		
		start = MPI_Wtime();
		// has even elements, which is guaranteed but last process
		
		// send to rank+1, last node do nothing
		if(rank!=size-1){
			MPI_Send(&array[former_alloc_num/2],former_alloc_num/2,MPI_INT,rank+1,0,MPI_COMM_WORLD);
			//MPI_Recv(&tmp1,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}
		// receive from rank-1, root node do nothing
		if(rank!=ROOT){
			MPI_Recv(temp_array,former_alloc_num/2,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);	  
		}
		// merge
		if(rank!=ROOT){
			for(int i=0,j=0,k=0;i<former_alloc_num;i++){
				if(j==former_alloc_num/2){
					sorted_array[i]=array[k];
					k++;
					continue;
				}
				else if(k==former_alloc_num/2){
					sorted_array[i]=temp_array[j];
					j++;
					continue;
				}
				if(temp_array[j]<array[k]){
					sorted_array[i]=temp_array[j];
					j++;
				}
				else{
					sorted_array[i]=array[k];
					k++;
				}
			}
		}
		// sorted array in sorted_array
		// original array in array
		// 0 to former_alloc_num/2 in array should be rear half part of sorted_array
		// first compare them 
		if(rank!=ROOT){
			if(memcmp(array, sorted_array+former_alloc_num/2, former_alloc_num*2)
				sorted = 0;
			memcpy(array, sorted_array+former_alloc_num/2, former_alloc_num*2);
		}
		if(rank!=ROOT){
			MPI_Send(sorted_array,former_alloc_num/2,MPI_INT,rank-1,0,MPI_COMM_WORLD);
		}
		// receive from rank+1, last node do nothing
		if(rank!=size-1){
			MPI_Recv(temp_array,former_alloc_num/2,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
		}
		if(rank!=size-1){
			if(memcmp(array+former_alloc_num/2, temp_array, former_alloc_num*2)
				sorted = 0;
			memcpy(array+former_alloc_num/2, temp_array, former_alloc_num*2);	
		}
		
    }
    //printf("rank: %2d %d %d\n", rank, array[0], array[1]);	
    MPI_Allreduce(&sorted,&sorted_temp,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    sorted = sorted_temp;
  }
#ifdef DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	//printall(array,alloc_num);
	MPI_Gatherv(array, alloc_num, MPI_INT, root_ptr, num_ptr, pos_ptr, MPI_INT, ROOT, MPI_COMM_WORLD);
	if(rank==ROOT){
		printall(root_ptr, N);
	}
#endif
	my_offset = rank*former_alloc_num*sizeof(int);
	start = MPI_Wtime();
	MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
	finish = MPI_Wtime();
	iotime += finish - start;
	MPI_Allreduce(&iotime,&io_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&commtime,&comm_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	io_all /= size;
	comm_all /= size;
	if(rank==ROOT)
		printf("iotime/t:%8.5lf/ncommtime/t:%8.5lf\n",io_all ,comm_all);
	MPI_Finalize();
	return 0;
}
