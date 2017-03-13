#pragma once
#include "types.h"

static void printSeg(struct InputData* seg) {
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

static void printComponents(struct MixtureComp* comps, int32_t nComp) {
  int i;
  printf("w:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].w);
  printf("\n");
  printf("mu:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].mu);
  printf("\n");
  printf("delta:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].delta);
  printf("\n");
  printf("sigmaSqF:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].sigmaSqF);
  printf("\n");
  printf("sigmaSqR:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].sigmaSqR);
  printf("\n");
  printf("se:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].se);
  printf("\n");
  printf("seF:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].seF);
  printf("\n");
  printf("seR:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].seR);
  printf("\n");
  printf("score:");
  for(i = 0; i < nComp; i++)
    printf(" %.4f", comps[i].score);
  printf("\n");
}

