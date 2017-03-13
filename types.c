#pragma once

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

static void free_input_data(struct InputData* data) {
  free(data->P.s);
  free(data->N.s);
  if(data->PC.s) {
    free(data->PC.s);
    free(data->NC.s);
  }
  if(data->U.s)
    free(data->U.s);
}

struct MixtureComp {
  double w, mu, delta, sigmaSqF, sigmaSqR, se, seF, seR, score;
};

struct MixtureResult {
  double BIC;
  int converge;
  int32_t nComp;
  struct MixtureComp* comps;
};

static struct MixtureResult* newMixtureResult(int32_t nComp) {
  struct MixtureResult* res = malloc(sizeof(struct MixtureResult));
  res->comps = malloc(sizeof(struct MixtureComp)*nComp);
  res->nComp = nComp;
  return res;
}

static void freeMixtureResult(struct MixtureResult* mr) {
  free(mr->comps);
  free(mr);
}

static struct MixtureResult* copyMixtureResult(struct MixtureResult* mr) {
  struct MixtureResult* copy = newMixtureResult(mr->nComp);
  copy->BIC = mr->BIC;
  copy->converge = mr->converge;
  int32_t i;
  for(i = 0; i < mr->nComp; i++)
     copy->comps[i] = mr->comps[i];
  return copy;
}

