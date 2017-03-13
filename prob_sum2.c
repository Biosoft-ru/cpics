static inline void FUNCTION_NAME (
  struct MixtureComp* comps, int32_t nComp,
  int32_t* Y, int32_t nY,
  struct SliceInterval* U, struct Interval* bounds,
  double* chi, double* eta) {

  memset(chi, 0, sizeof(double) * nComp);
  memset(eta, 0, sizeof(double) * nComp);

  for(int32_t i = 0; i < nY; i++) {
    double y = Y[i];
    double R[nComp];
    double Ruy[nComp];
    double sumR = 0;
    for(int32_t j = 0; j < nComp; j++) {
      double sigmaSq = SIGMASQ(comps[j]);
      double d = y - MU(comps[j]);
      d = d*d;
      double a = 1./(4*sigmaSq + d);
      double b = sigmaSq*a;
      double Rj = comps[j].w*b*b*sqrt(a);
      sumR += Rj;
      R[j] = Rj;
      Ruy[j] = 5*b*d*Rj;
    }
    sumR = 1./sumR;
    for(int32_t j = 0; j < nComp; j++) {
      chi[j] += R[j]*sumR;
      eta[j] += Ruy[j]*sumR;
    }
  }

  if(U->s == U->e)
    return;

  double P0 = 1;
  double H3[nComp], H2[nComp];
  memset(H3, 0, sizeof H3);
  memset(H2, 0, sizeof H2);
  for(int32_t j = 0; j < nComp; j++) {
    for(struct Interval* u = U->s; u < U->e; u++) {
      double sigma = sqrt(SIGMASQ(comps[j]));
      int32_t a = max(u->start, bounds->start);
      int32_t b = min(u->end, bounds->end) - 1;
      double aNorm = (a - MU(comps[j])) / sigma;
      double bNorm = (b - MU(comps[j])) / sigma;
      
      double cdfDiff = gsl_cdf_tdist_P(bNorm, 4) - gsl_cdf_tdist_P(aNorm, 4);
      H3[j] += cdfDiff;

      P0 -= comps[j].w * cdfDiff;
      
      double aSin = sin(atan(aNorm/2));
      double bSin = sin(atan(bNorm/2));

      H2[j] += cst*(fun2(bSin)-fun2(aSin));
    }
  }
  double phi = nY / P0;

  for(int32_t j = 0; j < nComp; j++) {
    double phiW = phi*comps[j].w;
    chi[j] += phiW * H3[j]; 
    eta[j] += phiW * SIGMASQ(comps[j]) * 4 * H2[j];
  }

}
