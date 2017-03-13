struct FunParams {
  int32_t nComp;
  double* w;
  double* mu;
  double* sigma;
};

static double plogp(double y, void* paramsP){
  struct FunParams* params = paramsP;
  double tmp = 0;
  for (int32_t j = 0; j < params->nComp; j++) {
    double yNorm = (y - params->mu[j])/params->sigma[j];
    tmp += params->w[j]* tdist4(yNorm) / params->sigma[j];      
  }
  return tmp*log(tmp);
}

static double intplogp(int32_t start, int32_t end, struct FunParams* params)
{
  gsl_integration_workspace * wk = gsl_integration_workspace_alloc(1000);
  double result, error;

  gsl_function FFF;
  FFF.function = &plogp;
  FFF.params = params;

  gsl_integration_qags (&FFF, start, end, 0, 1e-7, 1000, wk, &result, &error); 

  gsl_integration_workspace_free (wk);

  return result;
}

static double BIC_p_term(
  int32_t* Y, int32_t nY,
  struct FunParams* params,
  struct SliceInterval* U, struct Interval* bounds) {

  double bic = 0;
  for(int32_t i = 0; i < nY; i++) {
    double p = 0;
    for(int32_t j = 0; j < params->nComp; j++) {
      double yNorm = (Y[i] - params->mu[j])/params->sigma[j];
      p += params->w[j] * tdist4(yNorm)/params->sigma[j];
    }
    bic += log(p);
  }

  double P0 = 1;
  double integrate = 0;
  for(struct Interval* u = U->s; u < U->e; u++) {
    int32_t a = max(u->start, bounds->start);
    int32_t b = min(u->end, bounds->end) - 1;
    integrate += intplogp (a, b, params);
    for(int32_t j = 0; j < params->nComp; j++) {
      double aNorm = (a - params->mu[j]) / params->sigma[j];
      double bNorm = (b - params->mu[j]) / params->sigma[j];
      double cdfDiff = gsl_cdf_tdist_P(bNorm, 4) - gsl_cdf_tdist_P(aNorm, 4);
      P0 -= params->w[j] * cdfDiff;
    }
  }
  bic += nY/P0*integrate;

  return bic;
}

static double BIC(struct MixtureComp* comps, int nComp, struct InputData* seg) {
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;
  int32_t N = NF + NR;

  double bic = (1 - 5*nComp)*log(N)/2.0;

  double mu[nComp], sigma[nComp], w[nComp];
  struct FunParams params = {.nComp = nComp, .w = w, .mu = mu, .sigma = sigma};
  for(int32_t i = 0; i < nComp; i++) {
    mu[i] = comps[i].mu - 0.5*comps[i].delta;
    sigma[i] = sqrt(comps[i].sigmaSqF);
    w[i] = comps[i].w;
  }
  bic += BIC_p_term(seg->P.s, NF, &params, &seg->U, &seg->bounds);

  for(int32_t i = 0; i < nComp; i++) {
    mu[i] = comps[i].mu + 0.5*comps[i].delta;
    sigma[i] = sqrt(comps[i].sigmaSqR);
    w[i] = comps[i].w;
  }
  bic += BIC_p_term(seg->N.s, NR, &params, &seg->U, &seg->bounds);

  return bic;
}
