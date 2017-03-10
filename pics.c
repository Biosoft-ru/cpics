
#include <gsl/gsl_sf_gamma.h>
static double TDIST4_F;
static void init_tdist4() {
  TDIST4_F = gsl_sf_gamma(2.5)/gsl_sf_gamma(2)/sqrt(4*M_PI);
}
static inline double tdist4(double x) {
  double a = (1+0.25*x*x);
  a=a*a*sqrt(a);
  return TDIST4_F/a;
}

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

static inline void probMatrixF1(struct MixtureComp* comps, int32_t nComp, int32_t* Y, int32_t N,
  double* chi, double* nF, double* sF) {
  memset(nF, 0, sizeof(double) * nComp);
  memset(sF, 0, sizeof(double) * nComp);

  int i,j;
  double mu[nComp];
  for(j = 0; j < nComp; j++)
    mu[j] = comps[j].mu - 0.5*comps[j].delta;

  for(i = 0; i < N; i++) {
    double y = Y[i];
    double rF[nComp];
    double ruyF[nComp];
    double rFSum = 0;
    for(j = 0; j < nComp; j++) {
      double sigmaSq = comps[j].sigmaSqF;
      double a = y-mu[j];
      a = 4*sigmaSq + a*a;
      a = 1./a;
      double b = sigmaSq*a;
      double rFj = comps[j].w*b*b*sqrt(a);
      rFSum += rFj;
      rF[j] = rFj;
      ruyF[j] = 5*a*rFj;
    }
    rFSum = 1./rFSum;
    for(j = 0; j < nComp; j++) {
      chi[j] += rF[j]*rFSum;
      double ruyFj = ruyF[j]*rFSum;
      nF[j] += ruyFj;
      sF[j] += y * ruyFj;
    }
  }
}

static inline void probMatrixR1(struct MixtureComp* comps, int32_t nComp, int32_t* Y, int32_t N,
  double* chi, double* nR, double* sR) {
  memset(nR, 0, sizeof(double) * nComp);
  memset(sR, 0, sizeof(double) * nComp);

  int i,j;
  double mu[nComp];
  for(j = 0; j < nComp; j++)
    mu[j] = comps[j].mu + 0.5*comps[j].delta;

  for(i = 0; i < N; i++) {
    double y = Y[i];
    double rR[nComp];
    double ruyR[nComp];
    double rRSum = 0;
    for(j = 0; j < nComp; j++) {
      double sigmaSq = comps[j].sigmaSqR;
      double a = y-mu[j];
      a = 4*sigmaSq + a*a;
      a = 1./a;
      double b = sigmaSq*a;
      double rRj = comps[j].w*b*b*sqrt(a);
      rRSum += rRj;
      rR[j] = rRj;
      ruyR[j] = 5*a*rRj;
    }
    rRSum = 1./rRSum;
    for(j = 0; j < nComp; j++) {
      chi[j] += rR[j]*rRSum;
      double ruyRj = ruyR[j]*rRSum;
      nR[j] += ruyRj;
      sR[j] += y * ruyRj;
    }
  }
}

static void ECM1(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double chi[nComp];
  memset(chi, 0, sizeof chi);

  double nF[nComp], sF[nComp];
  probMatrixF1(comps, nComp, seg->P.s, NF, chi, nF, sF);

  double nR[nComp], sR[nComp];
  probMatrixR1(comps, nComp, seg->N.s, NR, chi, nR, sR);

  int j;
  for(j = 0; j < nComp; j++)
    comps[j].w = chi[j] / (NF + NR);

  int32_t dim = 2*nComp;
  double AA[dim*dim];
  memset(AA, 0, sizeof AA);
  double BB[dim];
  for(j = 0; j < nComp; j++) {
    AA[j*dim + j] = nF[j] + nR[j];
    AA[j*dim + (nComp+j)] = 0.5*(nR[j]-nF[j]);
    AA[(nComp+j)*dim + j] = nF[j]-nR[j];
    double ppp = 2*PriorParams.rho*(1.0/comps[j].sigmaSqF+1.0/comps[j].sigmaSqR);
    AA[(nComp+j)*dim + (nComp+j)] = -0.5*(nR[j]+nF[j]) - ppp;
    BB[j] = sF[j] + sR[j];
    BB[j+nComp] = sF[j] - sR[j] - ppp*PriorParams.xi;
  }

  /** Solve the system via LU decomposition **/
  gsl_matrix_view gsl_AA = gsl_matrix_view_array(AA, dim, dim);
  gsl_vector_view gsl_BB = gsl_vector_view_array(BB, dim);
  gsl_permutation *perm = gsl_permutation_alloc(dim);
  int sign;
  int status1 = gsl_linalg_LU_decomp(&gsl_AA.matrix, perm, &sign);
  double sol[dim];
  gsl_vector_view gsl_sol = gsl_vector_view_array(sol, dim);
  int status2 = gsl_linalg_LU_solve (&gsl_AA.matrix, perm, &gsl_BB.vector, &gsl_sol.vector);
  gsl_permutation_free(perm);

    
  /** Copy the new values **/
  /** Check that we could invert the matrix, otherwise we do not update the values **/
  if(status1==0 && status2==0)
    for(j=0;j<nComp;j++) {
      comps[j].mu = sol[j];
      comps[j].delta = sol[nComp + j];
    }
}

static inline void probMatrixF2(struct MixtureComp* comps, int32_t nComp, int32_t* Y, int32_t N,
  double* chiF, double* etaF) {
  memset(chiF, 0, sizeof(double) * nComp);
  memset(etaF, 0, sizeof(double) * nComp);

  int i,j;
  double mu[nComp];
  for(j = 0; j < nComp; j++)
    mu[j] = comps[j].mu - 0.5*comps[j].delta;

  for(i = 0; i < N; i++) {
    double y = Y[i];
    double rF[nComp];
    double ruyF[nComp];
    double rFSum = 0;
    for(j = 0; j < nComp; j++) {
      double sigmaSq = comps[j].sigmaSqF;
      double d = y - mu[j];
      d = d*d;
      double a = 1./(4*sigmaSq + d);
      double b = sigmaSq*a;
      double rFj = comps[j].w*b*b*sqrt(a);
      rFSum += rFj;
      rF[j] = rFj;
      ruyF[j] = 5*b*d*rFj;
    }
    rFSum = 1./rFSum;
    for(j = 0; j < nComp; j++) {
      chiF[j] += rF[j]*rFSum;
      etaF[j] += ruyF[j]*rFSum;
    }
  }
}

static inline void probMatrixR2(struct MixtureComp* comps, int32_t nComp, int32_t* Y, int32_t N,
  double* chiR, double* etaR) {
  memset(chiR, 0, sizeof(double) * nComp);
  memset(etaR, 0, sizeof(double) * nComp);

  int i,j;
  double mu[nComp];
  for(j = 0; j < nComp; j++)
    mu[j] = comps[j].mu + 0.5*comps[j].delta;

  for(i = 0; i < N; i++) {
    double y = Y[i];
    double rR[nComp];
    double ruyR[nComp];
    double rRSum = 0;
    for(j = 0; j < nComp; j++) {
      double sigmaSq = comps[j].sigmaSqR;
      double d = y - mu[j];
      d = d*d;
      double a = 1./(4*sigmaSq + d);
      double b = sigmaSq*a;
      double rRj = comps[j].w*b*b*sqrt(a);
      rRSum += rRj;
      rR[j] = rRj;
      ruyR[j] = 5*b*d*rRj;
    }
    rRSum = 1./rRSum;
    for(j = 0; j < nComp; j++) {
      chiR[j] += rR[j]*rRSum;
      etaR[j] += ruyR[j]*rRSum;
    }
  }
}

static void ECM2(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double chiF[nComp];
  double etaF[nComp];
  probMatrixF2(comps, nComp, seg->P.s, NF, chiF, etaF);

  double chiR[nComp];
  double etaR[nComp];
  probMatrixR2(comps, nComp, seg->N.s, NR, chiR, etaR);


  double cc = 2*PriorParams.alpha - 1;
  int j;
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

    ECM2(comps, nComp, seg);

    qsort(comps, nComp, sizeof(struct MixtureComp), &cmp_mixture_comp_mu);

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
      p += comps[j].w * tdist4(yNorm)/sigma;
    }
    bic += log(p);
  }
  for(i = 0; i < NR; i++) {
    double p = 0;
    for(j = 0; j < nComp; j++) {
      double sigma = sqrt(comps[j].sigmaSqR);
      double mu = comps[j].mu + comps[j].delta/2;
      double yNorm = (seg->N.s[i] - mu)/sigma;
      p += comps[j].w * tdist4(yNorm)/sigma;
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
    return 1; //Not enough reads
  
  initPara(res->comps, nComp, seg);
  int32_t iter = iterEM(res->comps, nComp, seg);
  res->converge = iter <= EMParams.B;
  
  //Handle the case of 'smallest weight too small'  
  int32_t i;
  for(i=0; i<nComp; i++)
    if (res->comps[i].w < 1.0/max(NF,NR))
      return 2; //too small weights
  
  res->BIC = BIC(res->comps, res->nComp, seg);

  return 0;
}  


static struct MixtureResult* fitPICS(struct InputData* seg) {
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
      if(best)
        freeMixtureResult(best);
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
