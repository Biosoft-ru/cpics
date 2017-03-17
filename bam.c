
struct BamReader {
  char* file_name;
  samFile* file;
  bam_hdr_t* header;
  hts_idx_t* idx;
};

static void open_bam_reader(char* file_name, struct BamReader* reader) {
  reader->file_name = strdup(file_name);
  reader->file = sam_open(file_name, "r");
  if(reader->file == 0) {
    fprintf(stderr, "Can not open bam file %s\n", file_name);
    exit(1);
  }
  reader->header = sam_hdr_read(reader->file);
  if(reader->header == 0) {
    fprintf(stderr, "Can not read header from bam file %s\n", file_name);
    exit(1);
  }
  reader->idx = sam_index_load(reader->file, reader->file_name);
  if(!reader->idx) {
    fprintf(stderr, "Can not find bam index file for %s\n", file_name);
    exit(1);
  }
}
static void close_bam_reader(struct BamReader* reader) {
  hts_idx_destroy(reader->idx);
  bam_hdr_destroy(reader->header);
  sam_close(reader->file);
  free(reader->file_name);
}

static int32_t count_mapped_reads(struct BamReader* reader, int32_t* chr_ids, int32_t n_chrs) {
  int32_t total_mapped = 0;
  for(int32_t i = 0; i < n_chrs; i++) {
    uint64_t mapped, unmapped;
    hts_idx_get_stat(reader->idx, chr_ids[i], &mapped, &unmapped);
    if( mapped > INT32_MAX - total_mapped) {
      fprintf(stderr, "Too many mapped reads in %s\n", reader->file_name);
      exit(1);
    }
    total_mapped += mapped;
  }
  return total_mapped;
}

/* Reads single chromosome data from bam file, the bam should be sorted and indexed
   chr_id = the chromosome to read
   P,N - arrays of 5' alignment ends for positive strand(P) and negative strand(N)
*/
static void read_chr_bam(struct BamReader* reader, int32_t chr_id,
   struct SliceInt32* P, struct SliceInt32* N) {
  
  size_t cP = 4;
  size_t nP = 0;
  P->s = malloc(sizeof(int32_t)*cP);
  if(P->s == 0)
    out_of_mem();

  size_t cN = 4;
  size_t nN = 0;
  N->s = malloc(sizeof(int32_t)*cN);
  if(N->s == 0)
    out_of_mem();

  hts_itr_t *iter = sam_itr_queryi(reader->idx, chr_id, 0, INT_MAX);
  if(!iter) {
    fprintf(stderr, "Can not navigate using index for %s\n", reader->file_name);
    exit(1);
  }

  bam1_t* rec = bam_init1();
  int result;
  while ((result = sam_itr_next(reader->file, iter, rec)) >= 0) {
    if(!(rec->core.flag & BAM_FUNMAP) && !(rec->core.qual < 10)) {
      if(bam_is_rev(rec)) {
        if(nN == cN) {
          N->s = realloc(N->s, sizeof(int32_t)*(cN *= 2));
          if(!N->s)
            out_of_mem();
        }
        N->s[nN++]= bam_endpos(rec); //one based 5' end
      } else {
        if(nP == cP) {
          P->s = realloc(P->s, sizeof(int32_t)*(cP *= 2));
          if(!P->s)
            out_of_mem();
        }
        P->s[nP++] = rec->core.pos + 1; //one based 5' end
      }
    }
  }

  if(result < -1) {
    fprintf(stderr, "Can not parse bam file %s\n", reader->file_name);
    exit(1);
  }

  bam_destroy1(rec);
  hts_itr_destroy(iter);

  P->e = P->s + nP;
  N->e = N->s + nN;
}
