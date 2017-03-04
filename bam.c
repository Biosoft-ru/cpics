
struct BamReader {
  samFile* file;
  bam_hdr_t* header;
  bam1_t* rec;
};

static int open_bam_reader(char* file_name, struct BamReader* reader) {
  reader->file = sam_open(file_name, "r");
  if(reader->file == 0)
    return 1;
  reader->header = sam_hdr_read(reader->file);
  if(reader->header == 0)
    return 2;
  reader->rec = 0;
  return 0;
}
static void close_bam_reader(struct BamReader* reader) {
  if(reader->rec != 0)
    bam_destroy1(reader->rec);
  sam_close(reader->file);
}

/* Reads single chromosome data from bam file, the bam should be sorted.
   chr_id = the chromosome to read
   P,N - arrays of 5' alignment ends for positive strand(P) and negative strand(N)
*/
static int read_chr_bam(struct BamReader* reader, int32_t chr_id,
   struct SliceInt32* P, struct SliceInt32* N) {
  
  size_t cP = 4;
  size_t nP = 0;
  P->s = malloc(sizeof(int32_t)*cP);
  if(P->s == 0)
    return 1;

  size_t cN = 4;
  size_t nN = 0;
  N->s = malloc(sizeof(int32_t)*cN);
  if(N->s == 0)
    return 1;


  if(reader->rec == 0) {
    reader->rec = bam_init1();
    if(sam_read1(reader->file, reader->header, reader->rec) < 0) {
      P->e = P->s;
      N->e = N->s;
      return 0;//EOF
    }
  }
  bam1_t* rec = reader->rec;

  int eof = 0;
  for(;;) {

    if(eof || rec->core.tid != chr_id) {
      P->e = P->s + nP;
      N->e = N->s + nN;
      return 0;
    }

    if(!(rec->core.flag & BAM_FUNMAP) && !(rec->core.qual < 10)) {
      if(bam_is_rev(rec)) {
        if(nN == cN)
          N->s = realloc(N->s, sizeof(int32_t)*(cN *= 2));
        N->s[nN++]= bam_endpos(rec); //one based 5' end
      } else {
        if(nP == cP)
          P->s = realloc(P->s, sizeof(int32_t)*(cP *= 2));
        P->s[nP++] = rec->core.pos + 1; //one based 5' end
      }
    }

    eof = sam_read1(reader->file, reader->header, rec) < 0;

  }
}
