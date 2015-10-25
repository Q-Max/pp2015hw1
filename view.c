#include <stdio.h>
#include <stdlib.h>
int main(int argc, char** argv){
	FILE *fp = fopen(argv[2],"rb");
	int count = atoi(argv[1]);
	int i;
	int *array = (int*)malloc(sizeof(int)*count);
	fread(array, sizeof(int), count, fp);
	for(i=0;i<count;i++)
		printf("%d ", array[i]);
	putchar('\n');
	return 0;
}
