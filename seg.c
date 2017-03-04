#include <stdio.h>
#include <stdlib.h>
#include <htslib/sam.h>

static int int32_cmp(const void* a, const void* b) {
  int32_t va = *(int32_t*)a;
  int32_t vb = *(int32_t*)b;
  return va - vb;
}

static inline int32_t min(int32_t a, int32_t b) { return a < b ? a : b; }
static inline int32_t max(int32_t a, int32_t b) { return a < b ? b : a; }


int step = 20;
int width = 250;
int min_reads = 2;
int min_reads_in_region = 3;
int min_l_region = 100;
char *exp_bam_file_name = 0, *ctrl_bam_file_name = 0, *unmappable_file_name = 0;

void usage() {
  fprintf(stderr, "Usage: cpics [-u <unmappable bed file>] [-c <control bam file>] <bam file> \n" );
}

int parse_args(int argc, char* argv[]) {
  int i;
  for(i = 1; i < argc; i++) {
    char* val = argv[i];
    if(*val == '-')  
      switch(*(val+1)) {
      case 'u':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -u option\n"); usage();
          return 1;
        }
        unmappable_file_name = argv[++i];
        break;
      case 'c':
        if(i+1>=argc) {
          fprintf(stderr, "No value for -c option\n"); usage();
          return 1;
        }
        ctrl_bam_file_name= argv[++i];
        break;
      default:
        fprintf(stderr, "Unknown option %s\n", val); usage();
        return 1;
      }
    else {
      if(exp_bam_file_name != 0) {
        fprintf(stderr, "extra argument %s\n", val); usage();
        return 1;
      }
      exp_bam_file_name = val;
    }
  }
  if(exp_bam_file_name == 0) {
    usage();
    return 1;
  }
  return 0;
}

struct SliceInt32 {
  int32_t *s, *e;// array of int32_t [s,e)
};

struct Interval {
  int32_t start, end; // [start,end) one based
};
struct SliceInterval {
  struct Interval *s, *e;// array of struct Interval [s,e)
};

struct InputData {
  char* chr;

  struct SliceInt32 P;//positive
  struct SliceInt32 N;//negative

  struct SliceInt32 PC;//positive control
  struct SliceInt32 NC;//negative control

  struct SliceInterval U;//intervals for unmappable regions
};

void free_input_data(struct InputData* data) {
  free(data->P.s);
  free(data->N.s);
  if(data->PC.s) {
    free(data->PC.s);
    free(data->NC.s);
  }
  if(data->U.s)
    free(data->U.s);
}

void processSeg(struct InputData* seg) {

  if(seg->P.e - seg->P.s < min_reads_in_region
  || seg->N.e - seg->N.s < min_reads_in_region)
    return;

  int32_t regionLen = *(seg->N.e-1) - *(seg->P.s);
  if(regionLen < min_l_region)
    return;

  int32_t* a;
  printf("%s \n", seg->chr);

  printf("F: ");
  for(a = seg->P.s; a < seg->P.e; a++)
    printf("%i ", *a);
  printf("\n");

  printf("R: ");
  for(a = seg->N.s; a < seg->N.e; a++)
    printf("%i ", *a);
  printf("\n");

  if(seg->PC.s) {
    printf("CF: ");
    for(a = seg->PC.s; a < seg->PC.e; a++)
      printf("%i ", *a);
    printf("\n");
    
    printf("CR: ");
    for(a = seg->NC.s; a < seg->NC.e; a++)
      printf("%i ", *a);
    printf("\n");
  }

  if(seg->U.s) {
    struct Interval* i;
    for(i = seg->U.s; i < seg->U.e; i++)
      printf("M: %i %i \n", i->start, i->end);
  }

}

void moveSliceInt32(struct SliceInt32* parent, struct SliceInt32* slice,
   int32_t minLoc, int32_t maxLoc) {
  while(slice->s < parent->e && *(slice->s) < minLoc)
    slice->s++;
  if(slice->e < slice->s)
    slice->e = slice->s;
  while(slice->e < parent->e && *(slice->e) <= maxLoc)
    slice->e++;
}

void prepareSeg(struct InputData* seg, struct InputData* data) {
  //the segment boundaries are defined by minimal coord on positive strand and maximal coord on negative strand
  int32_t minLoc = *seg->P.s;
  int32_t maxLoc = *(seg->N.e-1);

  int32_t* a = seg->P.s;
  while(a < data->P.e && *a <= maxLoc)
    ++a;
  seg->P.e = a;

  a = seg->N.e;
  while(a > data->N.s && *(a-1) >= minLoc)
    --a;
  seg->N.s = a;

  if(data->PC.s) {
    moveSliceInt32(&data->PC, &seg->PC, minLoc, maxLoc);
    moveSliceInt32(&data->NC, &seg->NC, minLoc, maxLoc);
  }

  if(data->U.s) {
    struct Interval* i = data->U.s;
    while(i < data->U.e && i->end <= minLoc)
      ++i;
    seg->U.s = i;

    if(seg->U.e < seg->U.s)
      seg->U.e = seg->U.s;

    i = seg->U.e;
    while(i < data->U.e && i->start <= maxLoc)
      ++i;
    seg->U.e = i;
  }
}

void process_chr(struct InputData* data) {
  int32_t* P = data->P.s;
  int32_t nP = data->P.e - data->P.s;
  int32_t* N = data->N.s;
  int32_t nN = data->N.e - data->N.s;
  if(nP == 0 || nN == 0)
    return;

  //P already sorted
  qsort(N, nN, sizeof(int32_t), int32_cmp);

  //scan chromosome with window [x-width,x+width] by step 
  //select windows with enough data and merge adjacent windows into segments

  //seg will hold all data in the current segment
  struct InputData seg = {0};
  seg.chr = data->chr;
  seg.PC.s = seg.PC.e = data->PC.s;
  seg.NC.s = seg.NC.e = data->NC.s;
  seg.U.s = seg.U.e = data->U.s;


  //sP will hold positive strand coordinates between x-width and x
  //sN will hold negative strand coordinates between x and x+width
  struct SliceInt32 sP, sN;
  sP.s = sP.e = data->P.s;
  sN.s = sN.e = data->N.s;
  int32_t segXEnd = -1;
  int32_t x;
  for(x = min(P[0],N[0]); x <= max(P[nP-1], N[nN-1]); x += step) {
    moveSliceInt32(&data->P, &sP, x-width, x);
    int32_t nPosL = sP.e - sP.s;

    moveSliceInt32(&data->N, &sN, x, x+width);
    int32_t nNegR = sN.e - sN.s;

    if(nPosL >= min_reads && nNegR >= min_reads) {
      if(segXEnd == -1) {
        seg.P.s = sP.s;
      } else if(x - segXEnd > 2*width) {
        prepareSeg(&seg, data);
        processSeg(&seg);
        seg.P.s = sP.s;
      }
      segXEnd = x;
      seg.N.e = sN.e;
    }

  }
  if(segXEnd != -1) {//Last segment
     prepareSeg(&seg, data);
     processSeg(&seg);
  }
}

struct BamReader {
  samFile* file;
  bam_hdr_t* header;
  bam1_t* rec;
};

int open_bam_reader(char* file_name, struct BamReader* reader) {
  reader->file = sam_open(file_name, "r");
  if(reader->file == 0)
    return 1;
  reader->header = sam_hdr_read(reader->file);
  if(reader->header == 0)
    return 2;
  reader->rec = 0;
  return 0;
}
void close_bam_reader(struct BamReader* reader) {
  if(reader->rec != 0)
    bam_destroy1(reader->rec);
  sam_close(reader->file);
}

/* Reads single chromosome data from bam file, the bam should be sorted.
   chr_id = the chromosome to read
   P,N - arrays of 5' alignment ends for positive strand(P) and negative strand(N)
*/
int read_chr(struct BamReader* reader, int32_t chr_id,
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

struct BedReader {
  FILE* file;
  struct Interval* last;
  char* lastChr;
};

int open_bed_reader(char* file_name, struct BedReader* reader) {
  reader->file = fopen(file_name, "r");
  if(!reader->file)
    return 1;
  reader->last = 0;
  reader->lastChr = 0;
  return 0;
}

void close_bed_reader(struct BedReader* reader) {
  fclose(reader->file);
  if(reader->last)
    free(reader->last);
  if(reader->lastChr)
    free(reader->lastChr);
}

int read_chr_bed(struct BedReader* reader, char* chr, struct SliceInterval* res) {
   int32_t capacity = 4;
   res->s = res->e = malloc(sizeof(struct Interval) * capacity);
   if(!reader->last) {
     reader->last = malloc(sizeof(struct BedReader));
     if(!reader->lastChr)
       reader->lastChr = malloc(32);
     int ret = fscanf(reader->file, "%31s %d %d", reader->lastChr, &(reader->last->start), &(reader->last->end));
     if(ret != 3 )
       return 0;
   }
   for(;;) {
     if(strcmp(reader->lastChr, chr) != 0)
       return 0;
     int32_t size = res->e - res->s;
     if(size == capacity) {
       res->s = realloc(res->s, sizeof(struct Interval) * (capacity *= 2));
       res->e = res->s + size;
     }
     res->e->start = reader->last->start;
     res->e->end = reader->last->end; 
     res->e++;
      
     int ret = fscanf(reader->file, "%31s %d %d", reader->lastChr, &(reader->last->start), &(reader->last->end));
     if(ret != 3) {
       free(reader->last);
       free(reader->lastChr);
       reader->last = 0;
       reader->lastChr = 0;
       return 0;
     }
   }
}

int main(int argc, char* argv[]) {
  if(parse_args(argc, argv) != 0)
    return 1;

  struct BamReader exp_bam_reader;
  open_bam_reader(exp_bam_file_name, &exp_bam_reader);
  int nChr = exp_bam_reader.header->n_targets;

  struct BamReader ctrl_bam_reader;
  if(ctrl_bam_file_name != 0) {
    open_bam_reader(ctrl_bam_file_name, &ctrl_bam_reader);
    if(ctrl_bam_reader.header->n_targets != nChr) {
      fprintf(stderr, "Distinct BAM headers in exp and ctrl files.\n");
      return 1;
    }
  }

  struct BedReader bed_reader;
  if(unmappable_file_name != 0) {
    if(open_bed_reader(unmappable_file_name, &bed_reader) != 0)
      return 1;
  }

  int i;
  for(i = 0; i < nChr; i++) {
    char* exp_chr = exp_bam_reader.header->target_name[i];
    char* ctrl_chr = ctrl_bam_reader.header->target_name[i];
    if(strcmp(exp_chr, ctrl_chr) != 0) {
      fprintf(stderr, "Distinct BAM headers in exp and ctrl files.\n");
      return 1;
    }

    struct InputData data = {0};
    data.chr = exp_chr;

    if(read_chr(&exp_bam_reader, i, &data.P, &data.N) != 0)
      return 1;

    if(ctrl_bam_file_name != 0) {
      if(read_chr(&ctrl_bam_reader, i, &data.PC, &data.NC) != 0)
        return 1;
    }

    if(unmappable_file_name != 0) {
      if(read_chr_bed(&bed_reader, exp_chr, &data.U) != 0)
        return 1;
    }

    process_chr(&data);
    free_input_data(&data);
  }

  close_bam_reader(&exp_bam_reader);
  if(ctrl_bam_file_name != 0)
    close_bam_reader(&ctrl_bam_reader);
  if(unmappable_file_name != 0)
    close_bed_reader(&bed_reader);

  return 0;
}
