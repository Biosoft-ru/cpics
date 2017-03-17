
char** load_chromosomes(int* n) {
  errno = 0;
  FILE* f = fopen(Opt.chr_file_name, "r");
  if(!f) {
    perror(Opt.chr_file_name);
    exit(1);
  }

  size_t capacity = 32;
  size_t size = 0;
  char** res = malloc(sizeof(char*)*capacity);
  if(!res)
    out_of_mem();

  char s[256];
  while(fgets(s, sizeof s, f)) {
    int l = 0;
    while(s[l] && s[l] != '\n' && s[l] != '\r')
      l++;
    if(l == 0) {
      fprintf(stderr, "Ignoring empty line in %s\n", Opt.chr_file_name);
      continue;
    }
    if(s[l]==0 && !feof(f)) {
      fprintf(stderr, "Chromosome name too long in file %s\n", Opt.chr_file_name);
      exit(1);
    }
    s[l] = 0;

    if(size == capacity) {
       capacity *= 2;
       res = realloc(res, sizeof(char*) * capacity);
       if(!res)
         out_of_mem();
    }
    res[size++] = strdup(s);
  }
  if(!feof(f)) {
    fprintf(stderr, "Can not read file %s\n", Opt.chr_file_name);
  }
  fclose(f);
  *n = size;
  return res;
}

void free_chromosomes(char** chrs, int n) {
  for(int i = 0; i < n; i++)
    free(chrs[i]);
  free(chrs);
}

static void match_chrs(struct BamReader* reader,  char* chrs[], int32_t n_chrs, int32_t* chr_ids) {
  for(int32_t i = 0; i < n_chrs; i++) {
    int chr_found = 0;
    for(int32_t j = 0; j < reader->header->n_targets; j++)
      if(strcmp(chrs[i], reader->header->target_name[j])==0) {
        chr_ids[i] = j;
        chr_found = 1;
        break;
      }
    if(!chr_found) {
      fprintf(stderr, "Can not find chromosome %s in bam file %s\n", chrs[i], reader->file_name);
      exit(1);
    }
  }
}
