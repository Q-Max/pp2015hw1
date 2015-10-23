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
			MPI_Finalize();
			exit(0);
		}
	}
	else if((alloc_num)%2){
		alloc_num--;
		MPI_File_seek(fp,(MPI_Offset)alloc_num*rank, MPI_SEEK_SET);
		if(rank==size-1){
			alloc_num = N - alloc_num * rank;
		}
		array = (int*)malloc(sizeof(int)*alloc_num);
		start = MPI_Wtime();
		MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
		finish = MPI_Wtime();
		printf("rank %2d io time: %lf\n", rank, finish-start);
	}
	else{
		// last process need deal with more inputs
    MPI_File_seek(fp, (MPI_Offset)rank*alloc_num*sizeof(int), MPI_SEEK_SET);
		if(rank==size-1){			
			alloc_num = N-rank*alloc_num;
		}	
		array = (int*)malloc(sizeof(int)*alloc_num);
		start = MPI_Wtime();
		MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
		finish = MPI_Wtime();
		printf("rank: %2d io time: %lf\n",rank,finish-start);	
  }
  // sheu
	// todo
	// sort and communicate with other
  //--------------------------------------------------------------------------------
  int tmp1,tmp2,sorted_temp; 
	int *tmp_ptr;
	int i,j,k,sorted=0;
	int nei_alloc;
	while(!sorted){
    sorted=1;
    for(i=0;i+1<alloc_num;i+=2){
      if(array[i]>array[i+1]){
        swap(&array[i],&array[i+1]);
        sorted = 0;
      }
    }
    printf("rank: %2d %d %d\n", rank, array[0], array[1]);	
    for(i=1;i<alloc_num;i+=2){
      if(i=alloc_num-1){
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
        /*if(array[i]>tmp2){
          swap(&array[i],&tmp2);
          sorted = 0;
        }*/
        /*
        // send to rank-1, last node do nothing
        if(rank!=size-1){
          MPI_Send(&array[i],1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
          //MPI_Recv(&tmp1,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        }
        // receive from rank+1, root node do nothing
        if(rank!=ROOT){
          MPI_Recv(&tmp1,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
        }
        if(array[i]>tmp1){
          swap(&array[i],&array[i+1]);
          sorted = 0;
        }
        // receive part, root node do nothing
        if(rank!=ROOT){
          MPI_Recv(&tmp1,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
          MPI_Send(&array[i],1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
        }
        */
      }
      else if(array[i]>array[i+1]){
        swap(&array[i],&array[i+1]);
        sorted = 0;
      }
    }
    printf("rank: %2d %d %d\n", rank, array[0], array[1]);
    MPI_Allreduce(&sorted,&sorted_temp,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
    sorted = sorted_temp;
  }

	MPI_Barrier(MPI_COMM_WORLD);
	printall(array,alloc_num);
	if(rank==ROOT){
		tmp_ptr = (int*)calloc(0,sizeof(int)*N);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gather(array,alloc_num,MPI_INT, tmp_ptr, alloc_num, MPI_INT, ROOT, MPI_COMM_WORLD);
	if(rank==ROOT){
		printall(tmp_ptr, N);
		free(tmp_ptr);
	}
	MPI_Finalize();
	exit(0);
  

	/*if(rank>=N){
		// sheu need to change!
		MPI_Finalize();
		exit(0);
	}
	else{
		if(N < size){
			size -= size - N;
		}
		alloc_num = N / size;
		if(N%size==0){
			//todo
			if(alloc_num==1){
				// all go to  last process
				if(rank!=size-1){
					MPI_Finalize();
					exit(0);
				}
				else{
					// in this situation, last node handle all input using singleoesort
					alloc_num = N;
					MPI_File_seek(fp,(MPI_Offset)0, MPI_SEEK_SET);
					array = (int*)malloc(sizeof(int)*alloc_num);
					start = MPI_Wtime();
					MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
					finish = MPI_Wtime();
					printf("rank %2d io time: %lf\n", rank, finish - start);
					singleOESort(array, alloc_num);
					printall(array, alloc_num);
					MPI_Finalize();
					exit(0);
			}
			MPI_File_seek(fp, (MPI_Offset)rank*alloc_num, MPI_SEEK_SET);
			array = (int*)malloc(sizeof(int)*alloc_num);
			start = MPI_Wtime();
			MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
			finish = MPI_Wtime();
			printf("rank: %2d io time: %lf\n",rank,finish-start);
		}
		else{
			// last process need to deal with more inputs
			if(rank==size-1){
				if(alloc_num%2){
					alloc_num--;
					MPI_File_seek(fp, (MPI_Offset)rank*alloc_num*sizeof(int), MPI_SEEK_SET);
					alloc_num = N-rank*alloc_num+size-1;
					// in some situation, last process need to deal with all process
				}
				else {
					MPI_File_seek(fp, (MPI_Offset)rank*alloc_num*sizeof(int), MPI_SEEK_SET);
					alloc_num = N-rank*alloc_num;
				}
				array = (int*)malloc(sizeof(int)*alloc_num);
				start = MPI_Wtime();
				MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
				finish = MPI_Wtime();
				printf("rank: %2d io time: %lf\n",rank,finish-start);

			}
			else{
				if(alloc_num%2)
					alloc_num--;
				MPI_File_seek(fp, (MPI_Offset)rank*alloc_num*sizeof(int), MPI_SEEK_SET);
				array = (int*)malloc(sizeof(int)*alloc_num);
				start = MPI_Wtime();
				MPI_File_read(fp, array, alloc_num, MPI_INT, &status);
				finish = MPI_Wtime();
				printf("rank: %2d io time: %lf\n",rank,finish-start);

			}*/
			/*alloc_num+=1;
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
			}*/
			//test
			/*
			insertionsort(array,alloc_num);
			printall(array,alloc_num);
			MPI_Finalize();
			exit(0);
			*/
	//	}
	
	// sheu
	// todo
	// sort and communicate with other
    //--------------------------------------------------------------------------------
/*	int time=0,sorted=0,temp;
	int *temp_array = malloc(sizeof(int)*(alloc_num+size));
	int *sorted_array = malloc(sizeof(int)*alloc_num);
	int *tmp_ptr;
	int i,j,k;
	int nei_alloc;
	while(!sorted){
		sorted = 1;
		printall(array,alloc_num);
		if(time==0){
			// from odd rank send to even rank
			if(rank%2){
				// odd rank
				if(rank==size-1){
					// do nothing
				}
				else{
					MPI_Send(&alloc_num,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
					MPI_Send(array,alloc_num,MPI_INT,rank+1,0,MPI_COMM_WORLD);
					MPI_Recv(&nei_alloc,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Recv(temp_array,nei_alloc,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					
				}
			}
			else{
				if(rank==ROOT){
					// do nothing
				}
				else{
					MPI_Recv(&nei_alloc,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Recv(temp_array,nei_alloc,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Send(&alloc_num,1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
					MPI_Send(array,alloc_num,MPI_INT,rank-1,0,MPI_COMM_WORLD);
				}
			}
			// merge
			// odd has smaller part
			if(rank%2){
				if(rank!=size-1){
					for(i=0,j=0,k=0;i<alloc_num;i++){
						if(j==nei_alloc){
							sorted_array[i]=array[k];
							k++;
							continue;
						}
						else if(k==alloc_num){
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
			}
			else{
			// even has larger part
				if(rank!=ROOT){
					for(i=alloc_num-1,j=nei_alloc-1,k=alloc_num-1;i>=0;i--){
						if(j==-1){
							sorted_array[i]=array[k];
							k--;
							continue;
						}
						else if(k==-1){
							sorted_array[i]=temp_array[j];
							j--;
							continue;
						}
						if(temp_array[j]>array[k]){
							sorted_array[i]=temp_array[j];
							j--;
						}
						else{
							sorted_array[i]=array[k];
							k--;
						}
					}
				}
			}
			if(rank!=ROOT&&(rank!=size-1||rank%2!=1)){
				tmp_ptr=array;
				array=sorted_array;
				sorted_array=tmp_ptr;
			}
			if(0!=memcmp((const void*)array,(const void*)sorted_array,(size_t)sizeof(int)*alloc_num))
				sorted=0;
			time=1;
		}
		if(time==1){
			// from even rank send to odd rank
			if(!(rank%2)){
				// even rank
				if(rank==size-1){
					// do nothing
				}
				else{
					MPI_Send(&alloc_num,1,MPI_INT,rank+1,0,MPI_COMM_WORLD);
					MPI_Send(array,alloc_num,MPI_INT,rank+1,0,MPI_COMM_WORLD);
					MPI_Recv(&nei_alloc,1,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Recv(temp_array,nei_alloc,MPI_INT,rank+1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					
				}
			}
			else{
				{
					MPI_Recv(&nei_alloc,1,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Recv(temp_array,nei_alloc,MPI_INT,rank-1,MPI_ANY_TAG,MPI_COMM_WORLD,&status);
					MPI_Send(&alloc_num,1,MPI_INT,rank-1,0,MPI_COMM_WORLD);
					MPI_Send(array,alloc_num,MPI_INT,rank-1,0,MPI_COMM_WORLD);
				}
			}
			// merge
//			printf("rank: %2d array",rank)
			// even has smaller part
			if(!(rank%2)){
				if(rank!=size-1){
					for(i=0,j=0,k=0;i<alloc_num;i++){
						if(j==nei_alloc){
							sorted_array[i]=array[k];
							k++;
							continue;
						}
						else if(k==alloc_num){
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
			}
			else{
			// odd has larger part
				for(i=alloc_num-1,j=nei_alloc-1,k=alloc_num-1;i>=0;i--){
					if(j==-1){
						sorted_array[i]=array[k];
						k--;
						continue;
					}
					else if(k==-1){
						sorted_array[i]=temp_array[j];
						j--;
						continue;
					}
					if(temp_array[j]>array[k]){
						sorted_array[i]=temp_array[j];
						j--;
					}
					else{
						sorted_array[i]=array[k];
						k--;
					}
				}
			}
			if(rank!=size-1||rank%2!=1){
				tmp_ptr=array;
				array=sorted_array;
				sorted_array=tmp_ptr;
			}
			if(0!=memcmp((const void*)array,(const void*)sorted_array,(size_t)sizeof(int)*alloc_num))
				sorted=0;

			time=0;
		}
		MPI_Allreduce(&sorted,&temp,1,MPI_INT,MPI_LAND,MPI_COMM_WORLD);
		sorted = temp;
	}
	printf("rank: %2d finish %d\n", rank, array[0]);
	
	//
	free(temp_array);
	if(rank==ROOT)
		temp_array = (int*)calloc(0,sizeof(int)*N);
	
	MPI_Gather(array,alloc_num,MPI_INT, temp_array, N, MPI_INT, ROOT, MPI_COMM_WORLD);
	if(rank==ROOT)
		printall(temp_array, N);
	MPI_Finalize();

	exit(0);
*/










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
/*	int num_per_node = N / size;
	int *local_arr = malloc (num_per_node * sizeof(int));
	// Ref:
	// MPI_Scatter (send_buf, send_count, send_type, recv_buf, recv_count, recv_type, root, comm)
	MPI_Scatter    (root_arr, num_per_node, MPI_INT, local_arr, num_per_node, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (rank == ROOT) {
		free (root_arr);
	}

	// Part 3: Calculate the sum of local_arr
	int sum = 0;
	//int i;
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
*/
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();

	return 0;
}
