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
			fprintf(stderr, "Usage: %s N input_file", argv[0]);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		return 0;
	}

	int N = atoi(argv[1]), alloc_num, former_alloc_num=0, last_alloc_num=0;
	const char *inName = argv[2];
	const char *outName = argv[3];
	double start, finish;
	int *root_ptr; // for root process (which rank == 0) only
  
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
	alloc_num = N/size;
	if(size>N){
		// if N < size, root take over
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
			printf("rank %2d io time: %lf\n", rank, finish - start);
			singleOESort(array, alloc_num);
			printall(array, alloc_num);
      MPI_File fh;
      MPI_File_open(MPI_COMM_WORLD, outName, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
      MPI_Offset my_offset = 0;
      MPI_File_write_at(fh, 0, array, alloc_num, MPI_INT, &status);
			MPI_Finalize();
			exit(0);
		}
	}
	else if(alloc_num==1){
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
			printf("rank %2d io time: %lf\n", rank, finish - start);
			singleOESort(array, alloc_num);
			printall(array, alloc_num);
      MPI_File fh;
      MPI_File_open(MPI_COMM_WORLD, outName, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
      MPI_Offset my_offset = 0;
      MPI_File_write_at(fh, 0, array, alloc_num, MPI_INT, &status);
			MPI_Finalize();
			exit(0);
		}
	}
	else if((alloc_num)%2){
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
		printf("rank %2d io time: %lf\n", rank, finish-start);
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
		printf("rank: %2d io time: %lf\n",rank,finish-start);	
  }
  //MPI_File_close(fp);
  // sheu
	// todo
	// sort and communicate with other
  //--------------------------------------------------------------------------------
  int tmp1,tmp2,sorted_temp;
  int *num_ptr, *pos_ptr; 
  int i,sorted=0;
#ifdef DEBUG
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
	while(!sorted){
    sorted=1;
    for(i=0;i+1<alloc_num;i+=2){
      if(array[i]>array[i+1]){
        swap(&array[i],&array[i+1]);
        sorted = 0;
      }
    }
    //printf("rank: %2d %d %d\n", rank, array[0], array[1]);	
    for(i=1;i<alloc_num;i+=2){
      if(i==alloc_num-1){
        // has even elements, which is guaranteed but last process
        
        // send to rank+1, last node do nothing
        if(rank!=size-1){
          MPI_Send(&array[i],1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
          //MPI_Recv(&tmp1,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
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
        
        // send to rank-1, root node do nothing
        if(rank!=ROOT){
          MPI_Send(&tmp1,1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
        }
        // receive from rank+1, last node do nothing
        if(rank!=size-1){
          MPI_Recv(&tmp2,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          array[i] = tmp2;
        }
        if(array[i]>tmp2&&rank!=size-1){
          swap(&array[i],&tmp2);
          sorted = 0;
        }
        /*if(rank==size-1){
          puts("NOOOO");
        }*/
        continue;
      }
      else if(array[i]>array[i+1]){
        swap(&array[i],&array[i+1]);
        sorted = 0;      
      }
      if(rank==size-1&&alloc_num%2&&i==alloc_num-2){
        // has odd elements, only in last process
        MPI_Recv(&tmp1,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);        
        if(array[0]<tmp1){
          swap(&array[0],&tmp1);
          sorted = 0;
        }
        MPI_Send(&tmp1,1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
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
  MPI_File fh;
  MPI_File_open(MPI_COMM_WORLD, outName, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
  MPI_Offset my_offset = rank*former_alloc_num*sizeof(int);
  MPI_File_write_at(fh, my_offset, array, alloc_num, MPI_INT, &status);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
