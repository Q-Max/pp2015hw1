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
	double start, finish, iotime = 0, commtime = 0, io_all, comm_all, cpu_all, cputime = 0, cpustart, cpufinish;
	int *root_ptr; // for root process (which rank == 0) only
  
	// Part 1: Read file
	/* Note: You should deal with cases where (N < size) in Homework 1 */
	int rc, i, trend = NOTSORTED;
	int *array;
	MPI_File fp;
	MPI_File fh;
	MPI_Offset my_offset;
	MPI_Request req;
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
			printf("iotime   : %8.5lf\ncommtime : %8.5lf\n",iotime,commtime);
			
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
			cpustart = MPI_Wtime();
			singleOESort(array, alloc_num);
			cpufinish = MPI_Wtime();
			cputime = cpufinish - cpustart;
#ifdef DEBUG
			printall(array, alloc_num);      
#endif
			my_offset = 0;
			start = MPI_Wtime();
			MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			iotime += finish - start;
			printf("iotime   : %8.5lfs\ncommtime : %8.5lfs\ncputime  : %8.5fs\n",iotime,commtime,cputime);
			MPI_Finalize();
			exit(0);
		}
	}
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
		array = (int*)malloc(sizeof(int)*alloc_num);
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
		array = (int*)malloc(sizeof(int)*alloc_num);
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
	cpustart = MPI_Wtime();
	while(!sorted){
		sorted=1;
		for(i=0;i+1<alloc_num;i+=2){
			if(array[i]>array[i+1]){
				swap(&array[i],&array[i+1]);
				sorted = 0;
			}
		}
		for(i=1;i<alloc_num;i+=2){
			if(i==alloc_num-1){
				start = MPI_Wtime();
				// has even elements, which is guaranteed but last process
			
				// send to rank+1, last node do nothing
				if(rank!=size-1){
					MPI_Isend(&array[i],1,MPI_INT,rank+1,0,MPI_COMM_WORLD,&req);
				}
				// receive from rank-1, root node do nothing
				if(rank!=ROOT){
					MPI_Recv(&tmp1,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);  
				}
				if(array[0]<tmp1&&rank!=ROOT){
					swap(&array[0],&tmp1);
					sorted = 0;
				}
				// tmp1 will always be smaller
				if(rank!=size-1)
					MPI_Wait(req, MPI_STATUS_IGNORE);
				// send to rank-1, root node do nothing
				if(rank!=ROOT){
					MPI_Send(&tmp1,1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
				}
				// receive from rank+1, last node do nothing
				if(rank!=size-1){
					MPI_Recv(&tmp2,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
				}
				if(array[i]>tmp2&&rank!=size-1){
					swap(&array[i],&tmp2);
					sorted = 0;
				}
				if(rank!=ROOT)
					MPI_Wait(req, MPI_STATUS_IGNORE);
				finish = MPI_Wtime();
				commtime += finish - start;
				continue;
			}
			else if(array[i]>array[i+1]){
				swap(&array[i],&array[i+1]);
				sorted = 0;
			}
			if(rank==size-1&&alloc_num%2&&i==alloc_num-2){
				start = MPI_Wtime();
				// has odd elements, only in last process
				MPI_Recv(&tmp1,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);        
				if(array[0]<tmp1){
					swap(&array[0],&tmp1);
					sorted = 0;
				}
				MPI_Isend(&tmp1,1,MPI_INT,rank-1,0,MPI_COMM_WORLD,&req);
				finish = MPI_Wtime();
				commtime += finish - start;
			}
		}
		start = MPI_Wtime();
		MPI_Allreduce(&sorted,&sorted_temp,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
		finish = MPI_Wtime();
		commtime += finish - start;
		sorted = sorted_temp;
	}
	cpufinish = MPI_Wtime();
	cputime = cpufinish - cpustart - commtime;
#ifdef DEBUG
	MPI_Barrier(MPI_COMM_WORLD);
	printall(array,alloc_num);
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
	MPI_Allreduce(&cputime,&cpu_all,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	io_all /= size;
	comm_all /= size;
	if(rank==ROOT)
		printf("iotime   : %8.5lfs\ncommtime : %8.5lfs\ncputime  : %8.5lfs(sum)\n",io_all ,comm_all, cpu_all);
	MPI_Finalize();
	return 0;
}
