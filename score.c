void computeScores(struct MixtureComp* comps, int32_t nComp, struct InputData* seg, int32_t N, int32_t Nc) {

  const double calpha=1.5;//percentile cutoff

  int32_t k;
  for(k=0;k<nComp;k++) {
    struct MixtureComp* c = comps + k;
    int32_t* p;

    double sigmaF = sqrt(c->sigmaSqF);
    double muF = c->mu - c->delta/2.;
    int32_t sumF = 0;
    for(p=seg->P.s;p<seg->P.e;p++)
      if(abs(*p - muF)/sigmaF < calpha)
        sumF++;

    double sigmaR = sqrt(c->sigmaSqR);
    double muR = c->mu + c->delta/2.;
    int32_t sumR = 0;
    for(p=seg->N.s;p<seg->N.e;p++)
      if(abs(*p - muR)/sigmaR < calpha)
        sumR++;

    double scoreF, scoreR;

    /** If we have control data, we can normalize the scores **/
    if(Nc > 0) { 
      int32_t sumCF = 0;
      for(p=seg->PC.s;p<seg->PC.e;p++)
        if(abs(*p - muF)/sigmaF < calpha)
          sumCF++;
      if(sumCF < 2)
        sumCF = 2;

      int32_t sumCR=0;
      for(p=seg->NC.s;p<seg->NC.e;p++)
        if(abs(*p - muR)/sigmaR < calpha)
          sumCR++;
      if(sumCR < 2)
        sumCR = 2;

      scoreF = sumF / (N*((double)sumCF)/Nc);
      scoreR = sumR / (N*((double)sumCR)/Nc);
    }
    else {
      scoreF = sumF;
      scoreR = sumR;
    }

    c->score = scoreF + scoreR;
  }
}
