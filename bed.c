
struct BedReader {
  FILE* file;
  struct Interval* last;
  char* lastChr;
};

static int open_bed_reader(char* file_name, struct BedReader* reader) {
  reader->file = fopen(file_name, "r");
  if(!reader->file)
    return 1;
  reader->last = 0;
  reader->lastChr = 0;
  return 0;
}

static void close_bed_reader(struct BedReader* reader) {
  fclose(reader->file);
  if(reader->last)
    free(reader->last);
  if(reader->lastChr)
    free(reader->lastChr);
}

static int read_chr_bed(struct BedReader* reader, char* chr, struct SliceInterval* res) {
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
      
     int32_t start, end;
     int ret = fscanf(reader->file, "%31s %d %d", reader->lastChr, &start, &end);
     reader->last->start = start + 1;
     reader->last->end = end + 1;
     if(ret != 3) {
       free(reader->last);
       free(reader->lastChr);
       reader->last = 0;
       reader->lastChr = 0;
       return 0;
     }
   }
}
