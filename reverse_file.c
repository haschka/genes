#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

int main(int argc, char** argv) {

  int i;

  int fd = open(argv[1], O_RDONLY);

  off_t size = lseek(fd, 0, SEEK_END);
  lseek(fd,0, SEEK_SET);

  char * buffer = (char*)malloc(sizeof(char)*size);

  read(fd, buffer, size);

  for(i = size-1;i>-1;i--) {
    putc(buffer[i],stdout);
  }
  return(0);
}
  
