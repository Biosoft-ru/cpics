
static int int32_cmp(const void* a, const void* b) {
  int32_t va = *(int32_t*)a;
  int32_t vb = *(int32_t*)b;
  return va - vb;
}

static int double_cmp(const void* a, const void* b) {
  double va = *(double*)a;
  double vb = *(double*)b;
  if(va < vb)
    return -1;
  if(va > vb)
    return 1;
  return 0;
}

static inline int32_t min(int32_t a, int32_t b) { return a < b ? a : b; }
static inline int32_t max(int32_t a, int32_t b) { return a < b ? b : a; }
