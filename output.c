void output(struct MixtureComp* comps, int32_t nComp, char* chr) {
  int i;
  for(i = 0;i < nComp; i++) {
    struct MixtureComp c = comps[i];
    if(c.score <= 10) continue;
    if(c.delta <= 50 || c.delta >= 300) continue;
    if(c.se >= 50 || c.se <= 0) continue; //!!! ==0 ???
    if(c.sigmaSqF > 22500 || c.sigmaSqR > 22500) continue;
    int32_t from = (int32_t) (c.mu - c.delta/2 - 3*c.seF);
    int32_t to = (int32_t) (c.mu + c.delta/2 + 3*c.seR);
    printf("%s\t%i\t%i\t.\t%f\n", chr, from - 1, to, c.score);
  }
}
