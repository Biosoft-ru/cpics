
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
  struct Interval bounds;//window bounds during segmentation
};

#include "pics.c"
#include "infmat.c"
#include "merge.c"
#include "score.c"
#include "output.c"

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

void swap_exp_ctrl(struct InputData* data) {
  struct SliceInt32 tmp;

  tmp = data->P;
  data->P = data->PC;
  data->PC = tmp;

  tmp = data->N;
  data->N = data->NC;
  data->NC = tmp;
}

void printSeg(struct InputData* seg) {
  int32_t* a;
  printf("%s\n", seg->chr);

  printf("F:");
  for(a = seg->P.s; a < seg->P.e; a++)
    printf(" %i", *a);
  printf("\n");

  printf("R:");
  for(a = seg->N.s; a < seg->N.e; a++)
    printf(" %i", *a);
  printf("\n");

  if(seg->PC.s) {
    printf("CF:");
    for(a = seg->PC.s; a < seg->PC.e; a++)
      printf(" %i", *a);
    printf("\n");
    
    printf("CR:");
    for(a = seg->NC.s; a < seg->NC.e; a++)
      printf(" %i", *a);
    printf("\n");
  }

  if(seg->U.s) {
    int32_t minBound = seg->bounds.start;
    int32_t maxBound = seg->bounds.end;
    struct Interval* i;
    for(i = seg->U.s; i < seg->U.e; i++)
      printf("M: %i %i\n", max(i->start,minBound), min(i->end, maxBound)-1);//print as one based, both inclusive
  }

}

void processSeg(struct InputData* seg, struct InputData* data, int32_t exp_read_count, int32_t ctrl_read_count) {
  if(seg->P.e - seg->P.s < min_reads_in_region
  || seg->N.e - seg->N.s < min_reads_in_region)
    return;

  //printSeg(seg);

  int32_t regionLen = *(seg->N.e-1) - *(seg->P.s);
  if(regionLen < min_l_region)
    return;

  struct MixtureResult* mix = fitPICS(seg);
  if(mix != 0) {
    double infMat[(5*mix->nComp-1)*(5*mix->nComp-1)];
    int flag = getInfMat(seg, mix->comps, mix->nComp, infMat);
    if(flag == 0) {
      memset(infMat, 0, sizeof infMat);// !!! zero infMat always???
      mergePeaks(mix->comps, &mix->nComp, infMat);
      computeScores(mix->comps, mix->nComp, seg, exp_read_count, ctrl_read_count);
      output(mix->comps, mix->nComp, seg->chr);
    }
    freeMixtureResult(mix);
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

void prepareSeg(struct InputData* seg, struct InputData* data, int32_t minLoc, int32_t maxLoc) {
  seg->bounds.start = minLoc;
  seg->bounds.end = maxLoc + 1;
  if(data->U.s) {
    struct Interval* i = data->U.s;
    while(i < data->U.e && i->end <= minLoc)
      ++i;
    seg->U.s = i;

    if(seg->U.e < seg->U.s)
      seg->U.e = seg->U.s;

    i = seg->U.e;
    while(i < data->U.e && i->start < maxLoc) // !!! incorrect, should be i->start <= maxLoc, but this way it is done in original pics 
      ++i;
    seg->U.e = i;
  }

  //the segment boundaries are defined by minimal coord on positive strand and maximal coord on negative strand
  maxLoc = *(seg->N.e-1);

  if(seg->U.e > seg->U.s && seg->U.s->start > minLoc && seg->U.s->start < *seg->P.s) {
    // !!! this is for compatibility with original pics only
    // !!! very odd part, in the original pics this part even can segfault
    minLoc = seg->U.s->start;
  } else {
    minLoc = *seg->P.s;
  }

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
}

void process_chr(struct InputData* data, int32_t exp_read_count, int32_t ctrl_read_count) {
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
  int32_t segXStart = -1;
  int32_t segXEnd = -1;
  int32_t x;
  for(x = min(P[0],N[0]); x <= max(P[nP-1], N[nN-1]); x += step) {
    moveSliceInt32(&data->P, &sP, x-width, x);
    int32_t nPosL = sP.e - sP.s;

    moveSliceInt32(&data->N, &sN, x, x+width);
    int32_t nNegR = sN.e - sN.s;

    if(nPosL >= min_reads && nNegR >= min_reads) {
      if(segXEnd == -1) { //no previous segment, start the first segment
        seg.P.s = sP.s;
        segXStart = x;
      } else if(x - segXEnd > 2*width) {//prev segment too far from current, report prev segment and start new one
        prepareSeg(&seg, data, segXStart - width, segXEnd + width);
        processSeg(&seg, data, exp_read_count, ctrl_read_count);
        seg.P.s = sP.s;
        segXStart = x;
      }
      segXEnd = x;
      seg.N.e = sN.e;
    }

  }
  if(segXEnd != -1) {//Last segment
     prepareSeg(&seg, data, segXStart - width, segXEnd + width);
     processSeg(&seg, data, exp_read_count, ctrl_read_count);
  }
}
