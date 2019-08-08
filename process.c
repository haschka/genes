#include<stdio.h>
#include<stdlib.h>
#include<string.h>

typedef struct {
  char** chromosomes;
  size_t* chromosome_size;
  int n_chromosomes;
  char** chromosome_name;
} gene;


gene read_gene_from_file(FILE* f, int n_chrome) {
  char* line = NULL;
  size_t n = 0;

  size_t name_size;
  char** chromosomes = (char**)malloc(sizeof(char*)*n_chrome);
  char** chromosome_name = (char**)malloc(sizeof(char*)*n_chrome);
  
  size_t* chromosome_size = (size_t*)malloc(sizeof(size_t)*n_chrome);

  int current_chromosome = -1;

  size_t alloc_size;

  gene retval;
  
  rewind(f);
  
  while(-1 != getline(&line, &n, f)) {
    
    if(line[0] == '>') {
      current_chromosome++;
      alloc_size = 10000;
      name_size = strlen(line)-2;
      chromosome_size[current_chromosome] = 0;
      chromosomes[current_chromosome] =
	(char*)malloc(sizeof(char)*alloc_size);
      chromosome_name[current_chromosome] =
	(char*)malloc(sizeof(char)*(name_size+1));
      memcpy(chromosome_name[current_chromosome],line+1,name_size);
      chromosome_name[current_chromosome][name_size] = 0;
    } else {
      // length of the string;
      n = strlen(line);
      // remove newline character
      n = n - 1;
      
      if (alloc_size < chromosome_size[current_chromosome] + n) {
	alloc_size += 10000;
	chromosomes[current_chromosome] =
	  (char*)realloc(chromosomes[current_chromosome],
			 sizeof(char*)*alloc_size);
      }
      memcpy(chromosomes[current_chromosome]
	     +chromosome_size[current_chromosome], line, n);
      chromosome_size[current_chromosome] += n;
    }
    n = 0;
  }
  retval.chromosomes = chromosomes;
  retval.chromosome_size = chromosome_size;
  retval.chromosome_name = chromosome_name;
  retval.n_chromosomes = n_chrome;
  return(retval);
}

long get_chromosome_index_by_name(gene gen, char* name) {

  long i;
  for(i = 0; i<gen.n_chromosomes; i++) {
    if(!strcmp(gen.chromosome_name[i],name)) return(i);
  }
  return(-1);
}

void extract_sequence_to_file(gene gen, size_t chromosome_index,
			      size_t start_position,
			      size_t end_position, FILE* f) {

  size_t i;
  
  size_t sequence_delta = end_position - start_position;

  size_t start_index = start_position - 1;
  size_t end_index = end_position - 1;

  char * inital_position = gen.chromosomes[chromosome_index]+start_index;
  
  for(i=0;i<sequence_delta;i++) {
    putc(inital_position[i],f);
  }
}


void apply_SNPs_from_vcf_on_gene(gene* gen, FILE* vcf_snp, double min_score) {

  char* line = NULL;
  size_t n = 0;

  char chromosome_name_buffer[80];
  size_t position;
  char buffer[80];
  char reference;
  char alternative;
  double quality;

  size_t base_index;
  size_t chromosome_index;

  long bases_modified = 0;
  
  while(-1 != getline(&line, &n, vcf_snp)) {
    if(line[0] != '#') {
      if( 6 == sscanf(line,"%s %lu %c %c %c %lf",
		      chromosome_name_buffer,
		      &position,
		      buffer,
		      &reference,
		      &alternative,
		      &quality) ) {
	
	chromosome_index =
	  (size_t)get_chromosome_index_by_name(gen[0],
					       chromosome_name_buffer);

	base_index = position - 1;
	
	if(quality > min_score) {
	  if(gen->chromosomes[chromosome_index][base_index] == reference) {
	    gen->chromosomes[chromosome_index][base_index] = alternative;
	    bases_modified++;
	  }
	}
      }
    }
  }
  printf("Applied SNP vcf, modified %lu bases!\n", bases_modified);	 
}	     

int main(int argc, char** argv) {

  // load data into ram:
  int n_chrome;
  char** chromosomes;

  FILE* genome_file = fopen(argv[1], "r");

  sscanf(argv[2],"%d",&n_chrome);

  read_gene_from_file(genome_file, n_chrome);

  return(0);
}
