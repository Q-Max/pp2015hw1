#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <unistd.h>
#include <fcntl.h>
int cmp(const void* a, const void* b);
int main(int argc, char** argv){

  if(argc!=4)
    exit(1);
  int numberToRead;
  if((numberToRead=atoi(argv[1]))==0)
    exit(1);
  struct stat st;
  int fd = 0;
  fd = open(argv[2],O_RDONLY);
  fstat(fd, &st);
  puts("QQ");
  if(fd != -1){
    if(st.st_size/sizeof(int)<numberToRead){
      puts("N is bigger than testcase in input file, read to end");
      numberToRead = st.st_size/sizeof(int);
    }
  }
  else{
    puts("exit");
    exit(1);
  }
  close(fd);
  FILE *fp = fopen(argv[2],"rb");
  int numbers[numberToRead];
  puts("GG");
  fread(numbers,sizeof(int),numberToRead,fp);
  qsort((void*)numbers,numberToRead,sizeof(int),cmp);
  for(int i=0;i<numberToRead;i++)
    printf("%d\n",numbers[i]);

  return 0;
}

int cmp(const void* a, const void* b){
  if(*(int*)a > *(int*)b)return 1;
  else if(*(int*)a < *(int*)b)return -1;
  return 0; 
}
