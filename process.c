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
  if (line != NULL) free(line);
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

int process_protein_exon_list(gene* gen, FILE* exon_list, FILE* output_file) {

  int i;
  
  char chromosome_name_buffer[80];

  size_t chromosome_index;
  
  char* line = NULL;
  size_t n = 0;

  char buffer[80];

  size_t start_position, start_index;
  size_t end_position, end_index;
  
  size_t length;
  int phase;

  char* chromosome;
  char* initial_position;
  
  rewind(exon_list);
  if (-1 != getline(&line, &n, exon_list)) {
    sscanf(line,"%s",chromosome_name_buffer);
    chromosome_index =
      (size_t)get_chromosome_index_by_name(gen[0], chromosome_name_buffer);
    chromosome = gen->chromosomes[chromosome_index];
  } else {
    fprintf(stderr, "Warning protein name not in file\n");
    return -1; // no protein name in file	
  }
  
  while(-1 != getline(&line, &n, exon_list)) {

    if( 5 == sscanf(line,"%s %lu %lu %lu %d",
		    buffer, &start_position, &end_position, &length, &phase)) {

      if(!strcmp(buffer,"Exon")) {

	start_index = start_position - 1;
	end_index = end_position - 1;

	if ( (length-1) == (end_index-start_index) ) {

	  initial_position = chromosome+start_index;

	  for (i = 0; i < length; i++) {
	    putc(initial_position[i], output_file);
	  }

	  fprintf(output_file,"\n");

	} else {
	  fprintf(stderr,"Waring lengths are not correct\n");
	  return -3; // length are not correct	  
	}
	
      } else {
	fprintf(stderr,"Warining line does not start with exon\n");
	return -2; // lines do not start with exon
      }
    } else {
      fprintf(stderr,"Warning format is not correct\n");
      return -4; // format not correct
    }
    
  }
  if (line != NULL) free(line);
  return 0; // file has been processed correctly
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
  
  if(line != NULL)free(line);
  printf("Applied SNP vcf, modified %lu bases!\n", bases_modified);
}	     

int main(int argc, char** argv) {

  // load data into ram:
  int n_chrome;
  char** chromosomes;

  FILE* genome_file = fopen(argv[1], "r");
  FILE* snp_vcf = fopen(argv[2], "r");
  FILE* exon_list = fopen(argv[3], "r");
  FILE* exon_out = fopen(argv[4], "w");

  gene gen;
  
  sscanf(argv[5],"%d",&n_chrome);

  gen = read_gene_from_file(genome_file, n_chrome);
  fclose(genome_file);
  apply_SNPs_from_vcf_on_gene(&gen, snp_vcf, 0.);
  process_protein_exon_list(&gen, exon_list, exon_out);
  
  return(0);
} 
