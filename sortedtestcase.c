#include <stdio.h>
#include <stdlib.h>
#include <string.h>
int main(int  argc, char** argv){
	int i, j, k;
	FILE *fp;
	int count = atoi(argv[1]);
	int trend = atoi(argv[2]);
	int num = atoi(argv[3]);
	char *outName = argv[4];
	int *array = (int*)malloc(sizeof(int)*count);
	for(i=0;i<count;i++){
		if(trend==0)
			array[i] = num;
		if(trend==1){
			array[i] = num + i;
		}
		if(trend==-1){
			array[i] = num - i;
		}
	}
	fp = fopen(outName, "wb");
	fwrite(array, sizeof(int), count, fp);
	return 0;
}
