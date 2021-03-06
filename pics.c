#include "types.c"

#include <gsl/gsl_sf_gamma.h>
static double TDIST4_F;
static double cst;

static void init_tdist4() {
  TDIST4_F = gsl_sf_gamma(2.5)/gsl_sf_gamma(2)/sqrt(4*M_PI);
  cst = gsl_sf_gamma(3.5)/gsl_sf_gamma(3.0)/M_SQRTPI;
}
static inline double tdist4(double x) {
  double a = (1+0.25*x*x);
  a=a*a*sqrt(a);
  return TDIST4_F/a;
}


static double fun0(double x)
{
  return pow(x,5)/5-pow(x,3)*2/3+x;
}

static double fun1(double x, int nu)
{
  return pow(1+x*x/4, -2.5);
}

static double fun2(double x)
{
  return pow(x,3)/3-pow(x,5)/5;
}

//Built in parameters, can not be changed by user.
static struct EMParams {
  int32_t maxK;//Maximal number of mixture components
  double tol;
  int B;//Maximum number of EM iterations
} EMParams = {
  .maxK = 15,
  .tol=1e-4,
  .B=100,
};

static struct PriorParams {
  double xi;
  double rho;
  double alpha;
  double beta;
} PriorParams = {
  .xi = 200,
  .rho = 1,
  .alpha = 20,
  .beta = 40000,
};

//PICS needs at least 3 reads on each strand per mixture component.
//This is distinct from similar parameter in segmentation procedure.
static int minReadsPerPeak = 3;

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

#define FUNCTION_NAME prob_sum1_F
#define MU(c) ((c).mu - 0.5*(c).delta)
#define SIGMASQ(c) ((c).sigmaSqF)
#include "prob_sum1.c"
#undef FUNCTION_NAME
#define FUNCTION_NAME prob_sum2_F
#include "prob_sum2.c"
#undef FUNCTION_NAME
#undef MU
#undef SIGMASQ

#define FUNCTION_NAME prob_sum1_R
#define MU(c) ((c).mu + 0.5*(c).delta)
#define SIGMASQ(c) ((c).sigmaSqR)
#include "prob_sum1.c"
#undef FUNCTION_NAME
#define FUNCTION_NAME prob_sum2_R
#include "prob_sum2.c"
#undef FUNCTION_NAME
#undef MU
#undef SIGMASQ

static void ECM1(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double chi[nComp];
  memset(chi, 0, sizeof chi);

  double nF[nComp], sF[nComp];
  double phiF;
  prob_sum1_F(comps, nComp, seg->P.s, NF, &seg->U, &seg->bounds, chi, nF, sF, &phiF);

  double nR[nComp], sR[nComp];
  double phiR;
  prob_sum1_R(comps, nComp, seg->N.s, NR, &seg->U, &seg->bounds, chi, nR, sR, &phiR);

  int j;
  for(j = 0; j < nComp; j++)
    comps[j].w = chi[j] / (phiF+phiR);

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

static void ECM2(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  double chiF[nComp];
  double etaF[nComp];
  prob_sum2_F(comps, nComp, seg->P.s, NF, &seg->U, &seg->bounds, chiF, etaF);

  double chiR[nComp];
  double etaR[nComp];
  prob_sum2_R(comps, nComp, seg->N.s, NR, &seg->U, &seg->bounds, chiR, etaR);


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

#include "bic.c"

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
