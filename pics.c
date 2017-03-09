
static struct EMParams {
  int32_t maxK;//Maximal number of mixture components
  double tol;
  int B;//Maximum number of EM iterations
  int mergePeaks;
  int mapCorrect;
} EMParams = {
  .maxK = 15,
  .tol=1e-4,
  .B=100,
  .mergePeaks=1,
  .mapCorrect=1
};

static struct PriorParams {
  double xi;
  double rho;
  double alpha;
  double beta;
  double dMu;//for histones
  int PExi;
} PriorParams = {
  .xi = 200,
  .rho = 1,
  .alpha = 20,
  .beta = 40000,
  .dMu = 0,
  .PExi = 0
};

//??? distinct from similar parameters in segmentation procedure
static int minReadsPerPeak = 3;
static int minReadsPerRegion = 4;

struct MixtureComp {
  double w, mu, delta, sigmaSqF, sigmaSqR, se, seF, seR, score;
};

struct MixtureResult {
  double BIC;
  int converge;
  int32_t nComp;
  struct MixtureComp* comps;
};

void printComponents(struct MixtureComp* comps, int32_t nComp) {
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

static void initPara(struct MixtureComp* comps, int32_t nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double varF = gsl_stats_int_variance(seg->P.s,1,NF)/nComp;
  double varR = gsl_stats_int_variance(seg->N.s,1,NR)/nComp;  
  
  int32_t i;
  for (i=0; i<nComp; i++) {
    struct MixtureComp* mc = comps + i;
    double muF = gsl_stats_int_quantile_from_sorted_data(seg->P.s, 1, NF, (2*i+1.0)/(2*nComp));
    double muR = gsl_stats_int_quantile_from_sorted_data(seg->N.s, 1, NR, (2*i+1.0)/(2*nComp));
    mc->mu = (muF+muR)/2;
    mc->w = 1.0/nComp;
    mc->delta = 150;
    mc->sigmaSqF = varF;
    mc->sigmaSqR = varR;
  }
}

static void ECM1(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double chi[nComp];
  memset(chi, 0, sizeof chi);
  double nF[nComp], sF[nComp];
  memset(nF, 0, sizeof nF);
  memset(sF, 0, sizeof sF);
  //printf("rF:\n");
  //printf("ruyF:\n");
  int i,j;
  for(i = 0; i < NF; i++) {
    double rF[nComp];
    double ruyF[nComp];
    double rFSum = 0;
    for(j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double fMu = c->mu - c->delta/2;
      double sigma = sqrt(c->sigmaSqF);
      double yNormF = (seg->P.s[i] - fMu)/sigma;
      rF[j] = c->w * gsl_ran_tdist_pdf(yNormF, 4) /sigma;
      rFSum += rF[j];
      ruyF[j] = 5.0/(4 + yNormF*yNormF);
    }
    for(j = 0; j < nComp; j++) {
      rF[j] /= rFSum;
      //printf(" %.4f", rF[j]);
      ruyF[j] *= rF[j];
    }
    for(j = 0; j < nComp; j++) {
      chi[j] += rF[j];
      nF[j] += ruyF[j];
      //printf(" %.4f", ruyF[j]);
      sF[j] += seg->P.s[i] * ruyF[j];
    }
    //printf("\n");
  }
  for(j = 0; j < nComp; j++) {
      nF[j] /= comps[j].sigmaSqF;
      sF[j] /= comps[j].sigmaSqF;
  }

  double nR[nComp], sR[nComp];
  memset(nR, 0, sizeof nR);
  memset(sR, 0, sizeof sR);
  //printf("rR:\n");
  for(i = 0; i < NR; i++) {
    double rR[nComp];
    double ruyR[nComp];
    double rRSum = 0;
    for(j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double rMu = c->mu + c->delta/2;
      double sigma = sqrt(c->sigmaSqR);
      double yNormR = (seg->N.s[i] - rMu)/sigma;
      rR[j] = c->w * gsl_ran_tdist_pdf(yNormR, 4)/sigma;
      rRSum += rR[j];
      ruyR[j] = 5.0/(4 + yNormR*yNormR);
    }
    for(j = 0; j < nComp; j++) {
      rR[j] /= rRSum;
      //printf(" %.4f", rR[j]);
      ruyR[j] *= rR[j];
    }
    for(j = 0; j < nComp; j++) {
      chi[j] += rR[j];
      nR[j] += ruyR[j];
      sR[j] += seg->N.s[i] * ruyR[j];
    }
    //printf("\n");
  }
  for(j = 0; j < nComp; j++) {
      nR[j] /= comps[j].sigmaSqR;
      sR[j] /= comps[j].sigmaSqR;
  }

  for(j = 0; j < nComp; j++)
    comps[j].w = chi[j] / (NF + NR);

 //printf("nF:");
 //for(j = 0; j < nComp; j++)
 //  printf(" %.4f", nF[j]);
 //printf("\n");

 //printf("nR:");
 //for(j = 0; j < nComp; j++)
 //  printf(" %.4f", nR[j]);
 //printf("\n");

 //printf("sF:");
 //for(j = 0; j < nComp; j++)
 //  printf(" %.4f", sF[j]);
 //printf("\n");

 //printf("sR:");
 //for(j = 0; j < nComp; j++)
 //  printf(" %.4f", sR[j]);
 //printf("\n");

 //printf("ppp:");
  int32_t dim = 2*nComp;
  double AA[dim*dim];
  memset(AA, 0, sizeof AA);
  double BB[dim];
  for(j = 0; j < nComp; j++) {
    AA[j*dim + j] = nF[j] + nR[j];
    AA[j*dim + (nComp+j)] = 0.5*(nR[j]-nF[j]);
    AA[(nComp+j)*dim + j] = nF[j]-nR[j];
    double ppp = 2*PriorParams.rho*(1.0/comps[j].sigmaSqF+1.0/comps[j].sigmaSqR);
    //printf(" %.4f", ppp);
    AA[(nComp+j)*dim + (nComp+j)] = -0.5*(nR[j]+nF[j]) - ppp;
    BB[j] = sF[j] + sR[j];
    BB[j+nComp] = sF[j] - sR[j] - ppp*PriorParams.xi;
  }
  //printf("\n");

  /** Solve the system via LU decomposition **/
  gsl_matrix_view gsl_AA = gsl_matrix_view_array(AA, dim, dim);
  gsl_vector_view gsl_BB = gsl_vector_view_array(BB, dim);
  gsl_permutation *perm = gsl_permutation_alloc(dim);
  int sign;
  int status1 = gsl_linalg_LU_decomp(&gsl_AA.matrix, perm, &sign);
  double sol[dim];
  gsl_vector_view gsl_sol = gsl_vector_view_array(sol, dim);
  int status2 = gsl_linalg_LU_solve (&gsl_AA.matrix, perm, &gsl_BB.vector, &gsl_sol.vector);
    
  /** Copy the new values **/
  /** Check that we could invert the matrix, otherwise we do not update the values **/
  if(status1==0 && status2==0)
    for(j=0;j<nComp;j++) {
      comps[j].mu = sol[j];
      comps[j].delta = sol[nComp + j];
    }
}

static void ECM2(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double chiF[nComp];
  memset(chiF, 0, sizeof chiF);
  double etaF[nComp];
  memset(etaF, 0, sizeof etaF);
  int i,j;
  for(i = 0; i < NF; i++) {
    double rF[nComp];
    double ruyF[nComp];
    double rFSum = 0;
    for(j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double fMu = c->mu - c->delta/2;
      double sigma = sqrt(c->sigmaSqF);
      double yF = (seg->P.s[i] - fMu);
      double yNormF = yF/sigma;
      rF[j] = c->w * gsl_ran_tdist_pdf(yNormF, 4)/sigma;
      rFSum += rF[j];
    }
    for(j = 0; j < nComp; j++) {
      rF[j] /= rFSum;
      struct MixtureComp* c = comps + j;
      double fMu = c->mu - c->delta/2;
      double sigma = sqrt(c->sigmaSqF);
      double yF = (seg->P.s[i] - fMu);
      double yNormF = yF/sigma;
      ruyF[j] = (rF[j]*5.0/(4 + yNormF*yNormF))*yF*yF;
    }
    for(j = 0; j < nComp; j++) {
      chiF[j] += rF[j];
      etaF[j] += ruyF[j];
    }
  }

  double chiR[nComp];
  memset(chiR, 0, sizeof chiR);
  double etaR[nComp];
  memset(etaR, 0, sizeof etaR);
  for(i = 0; i < NR; i++) {
    double rR[nComp];
    double ruyR[nComp];
    double rRSum = 0;
    for(j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double rMu = c->mu + c->delta/2;
      double sigma = sqrt(c->sigmaSqR);
      double yNormR = (seg->N.s[i] - rMu)/sigma;
      rR[j] = c->w * gsl_ran_tdist_pdf(yNormR, 4)/sigma;
      rRSum += rR[j];
    }
    for(j = 0; j < nComp; j++) {
      rR[j] /= rRSum;
      struct MixtureComp* c = comps + j;
      double rMu = c->mu + c->delta/2;
      double sigma = sqrt(c->sigmaSqR);
      double yR = (seg->N.s[i] - rMu);
      double yNormR = yR/sigma;
      ruyR[j] = (rR[j]*5.0/(4 + yNormR*yNormR))*yR*yR;
    }
    for(j = 0; j < nComp; j++) {
      chiR[j] += rR[j];
      etaR[j] += ruyR[j];
    }
  }


  double cc = 2*PriorParams.alpha - 1;
  for(j=0;j<nComp;j++)
  {
    double dd = comps[j].delta - PriorParams.xi;
    dd = PriorParams.rho*dd*dd + 2*PriorParams.beta;
    
    comps[j].sigmaSqF = (etaF[j]+dd)/(cc+chiF[j]);
    comps[j].sigmaSqR = (etaR[j]+dd)/(cc+chiR[j]);
  }
}

static int cmp_mixture_comp_mu(const void* a, const void* b) {
  struct MixtureComp* ma = (struct MixtureComp*) a;
  struct MixtureComp* mb = (struct MixtureComp*) b;
  if(ma->mu < mb->mu)
    return -1;
  if(ma->mu > mb->mu)
    return 1;
  return 0;
}

static int32_t iterEM(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int iter;
  for(iter = 1; iter <= EMParams.B; iter++) {
    double old_mu[nComp];
    int i;
    for(i = 0; i < nComp; i++)
      old_mu[i] = comps[i].mu;

    ECM1(comps, nComp, seg);
    //printf("ECM1 ITER%i:\n", iter);
    //printComponents(comps, nComp);

    ECM2(comps, nComp, seg);

    qsort(comps, nComp, sizeof(struct MixtureComp), &cmp_mixture_comp_mu);

    //printf("ECM2 ITER%i:\n", iter);
    //printComponents(comps, nComp);

    double delta_sum = 0;
    for(i = 0; i < nComp; i++)
      delta_sum += fabs(comps[i].mu - old_mu[i]);   
    if(delta_sum < EMParams.tol)
      break;
  }
  return iter;
}

static double BIC(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;
  int32_t N = NF + NR;

  double bic = (1 - 5*nComp)*log(N)/2.0;
  int32_t i,j;
  for(i = 0; i < NF; i++) {
    double p = 0;
    for(j = 0; j < nComp; j++) {
      double sigma = sqrt(comps[j].sigmaSqF);
      double mu = comps[j].mu - comps[j].delta/2;
      double yNorm = (seg->P.s[i] - mu)/sigma;
      p += comps[j].w * gsl_ran_tdist_pdf(yNorm, 4)/sigma;// ??? /sigma
    }
    bic += log(p);
  }
  for(i = 0; i < NR; i++) {
    double p = 0;
    for(j = 0; j < nComp; j++) {
      double sigma = sqrt(comps[j].sigmaSqR);
      double mu = comps[j].mu + comps[j].delta/2;
      double yNorm = (seg->N.s[i] - mu)/sigma;
      p += comps[j].w * gsl_ran_tdist_pdf(yNorm, 4)/sigma;// ??? /sigma
    }
    bic += log(p);
  }

  return bic;
}

static int fitModel(int nComp, struct InputData* seg, struct MixtureResult* res) {

  res->nComp = nComp;

  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;
  
  if(nComp*minReadsPerPeak > min(NF, NR))
  {
    //printf("NOT ENOUGH READS\n");
    return 1; //Not enough reads
  }
  
  initPara(res->comps, nComp, seg);
  //printf("INIT:\n");
  //printComponents(res->comps, nComp);
  int32_t iter = iterEM(res->comps, nComp, seg);
  res->converge = iter <= EMParams.B;
  
  //Handle the case of 'smallest weight too small'  
  int32_t i;
  for(i=0; i<nComp; i++)
    if (res->comps[i].w < 1.0/max(NF,NR)) {
      //printf("TOO SMALL WEIGHTS\n");
      return 2; //too small weights
    }
  
  res->BIC = BIC(res->comps, res->nComp, seg);
  //printf("NCOMP:%i,BIC:%.4f\n", nComp, res->BIC);

  return 0;
}  


static struct MixtureResult* fitPICS(struct InputData* seg, int32_t total_reads_exp, int32_t total_reads_ctrl) {
  int32_t range = max(*(seg->P.e-1), *(seg->N.e-1)) - min(*(seg->P.s), *(seg->N.s));
  int32_t maxComp = (int32_t) 1 + range/(PriorParams.xi+4.0*sqrt(PriorParams.beta/PriorParams.alpha));
  if(maxComp > EMParams.maxK)
    maxComp = EMParams.maxK;

  struct MixtureResult* best = 0;
  double bestBIC = -INFINITY;
  int32_t nDecreaseBIC = 0;
  int32_t nComp;
  for(nComp = 1; nComp <= maxComp && nDecreaseBIC < 2; nComp++) {
    struct MixtureResult* res = newMixtureResult(nComp);
    if(fitModel(nComp, seg, res) == 0 && res->BIC > bestBIC) {
      best = copyMixtureResult(res);
      bestBIC = res->BIC;
      nDecreaseBIC = 0;
    }
    else
      nDecreaseBIC++;
    freeMixtureResult(res);
  }
  return best;
}
