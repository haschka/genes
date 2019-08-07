#include<stdio.h>
#include<stdlib.h>

typedef struct {
  char** chromosomes;
  size_t* chromosome_size;
  int n_chromosomes;
} gene;


gene read_gene_from_file(FILE* f, int n_chrome) {
  char* line = NULL;
  size_t n = 0;

  char** chromosomes = (char**)malloc(sizeof(char*)*n_chrome);

  size_t* chromosome_size = (size_t*)malloc(sizeof(size_t)*n_chrome);

  int current_chromosome = -1;

  size_t alloc_size;

  gene retval;
  
  rewind(f);
  
  while(-1 != getline(&line, &n, f)) {

    if(line[0] == '>') {
      current_chromosome++;
      alloc_size = 10000;
      chromosome_size[current_chromosome] = 0;
      chromosomes[current_chromosome] =
	(char*)malloc(sizeof(char*)*alloc_size);

    } else {
      if (alloc_size < chromosome_size + n) {
	alloc_size += 10000;
	chromosomes[current_chromosome] =
	  (char*)realloc(chromosomes[current_chromosome],
			 sizeof(char*)*alloc_size);
      }
      memcpy(chromosomes[current_chromosome]
	     +chromosome_size[current_chromosome], line, n);
      chromosome_size[current_chromosome] += n;
    }	
  }
  gene.chromosomes = chromosomes;
  gene.chromosome_size = chromosome_size;
  gene.n_chromosomes = n_chrome;
}
    

int main(int argc, char** argv) {

  // load data into ram:

  char** chromosomes = (char*)malloc(sizeof(char*)*n_chrome);

  FILE* genome_file = fopen(argv[1], "r");

  int n_chrome;

  sscanf(argv[2],"%d",n_chrome);

  read_gene_from_file(genome_file, n_chrome);
}
