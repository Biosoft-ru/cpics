#include "types.c"
#include "pics.c"
#include "infmat.c"
#include "merge.c"
#include "score.c"
#include "output.c"

static void swap_exp_ctrl(struct InputData* data) {
  struct SliceInt32 tmp;

  tmp = data->P;
  data->P = data->PC;
  data->PC = tmp;

  tmp = data->N;
  data->N = data->NC;
  data->NC = tmp;
}

static int32_t getSegmentLen(struct InputData* seg) {
  return *(seg->N.e-1) - *(seg->P.s);
}

static void processSeg(struct InputData* seg, struct InputData* data, int32_t exp_read_count, int32_t ctrl_read_count) {
  if(seg->P.e - seg->P.s < Opt.min_reads_in_region
  || seg->N.e - seg->N.s < Opt.min_reads_in_region)
    return;

  if(getSegmentLen(seg) < Opt.min_l_region)
    return;

  struct MixtureResult* mix = fitPICS(seg);
  if(mix != 0) {
    double infMat[(5*mix->nComp-1)*(5*mix->nComp-1)];
    int flag = getInfMat(seg, mix->comps, mix->nComp, infMat);
    if(flag == 0) {
      mergePeaks(mix->comps, &mix->nComp, infMat);
      computeScores(mix->comps, mix->nComp, seg, exp_read_count, ctrl_read_count);
      output(mix->comps, mix->nComp, seg->chr);
    }
    freeMixtureResult(mix);
  }
}

static void moveSliceInt32(struct SliceInt32* parent, struct SliceInt32* slice,
   int32_t minLoc, int32_t maxLoc) {
  while(slice->s < parent->e && *(slice->s) < minLoc)
    slice->s++;
  if(slice->e < slice->s)
    slice->e = slice->s;
  while(slice->e < parent->e && *(slice->e) <= maxLoc)
    slice->e++;
}

static void prepareSeg(struct InputData* seg, struct InputData* data, int32_t minLoc, int32_t maxLoc) {
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
    while(i < data->U.e && i->start <= maxLoc)
      ++i;
    seg->U.e = i;
  }

  //the segment boundaries are defined by minimal coord on positive strand and maximal coord on negative strand
  minLoc = *seg->P.s;
  maxLoc = *(seg->N.e-1);

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

static void process_chr(struct InputData* data, int32_t exp_read_count, int32_t ctrl_read_count) {
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
  for(x = min(P[0],N[0]); x <= max(P[nP-1], N[nN-1]); x += Opt.step) {
    moveSliceInt32(&data->P, &sP, x-Opt.width, x);
    int32_t nPosL = sP.e - sP.s;

    moveSliceInt32(&data->N, &sN, x, x+Opt.width);
    int32_t nNegR = sN.e - sN.s;

    if(nPosL >= Opt.min_reads_in_window && nNegR >= Opt.min_reads_in_window) {
      if(segXEnd == -1) { //no previous segment, start the first segment
        seg.P.s = sP.s;
        segXStart = x;
      } else if(x - segXEnd > 2*Opt.width) {//prev segment too far from current, report prev segment and start new one
        prepareSeg(&seg, data, segXStart - Opt.width, segXEnd + Opt.width);
        processSeg(&seg, data, exp_read_count, ctrl_read_count);
        seg.P.s = sP.s;
        segXStart = x;
      }
      segXEnd = x;
      seg.N.e = sN.e;
    }

  }
  if(segXEnd != -1) {//Last segment
     prepareSeg(&seg, data, segXStart - Opt.width, segXEnd + Opt.width);
     processSeg(&seg, data, exp_read_count, ctrl_read_count);
  }
}

static void process_chr_timed(struct InputData* data, int32_t exp_read_count, int32_t ctrl_read_count) {
  fprintf(stderr, "Calling peaks...");
  clock_t start_time = clock();
  process_chr(data, exp_read_count, ctrl_read_count);
  double run_time = (double)(clock() - start_time)/CLOCKS_PER_SEC;
  fprintf(stderr, "(%lfsec)\n", run_time);
}
