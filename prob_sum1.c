static inline void FUNCTION_NAME (
  struct MixtureComp* comps, int32_t nComp,
  int32_t* Y, int32_t nY,
  struct SliceInterval* U, struct Interval* bounds,
  double* chi, double* n, double* s, double* phi) {

  memset(n, 0, sizeof(double) * nComp);
  memset(s, 0, sizeof(double) * nComp);

  for(int32_t i = 0; i < nY; i++) {
    double y = Y[i];
    double R[nComp];
    double Ruy[nComp];
    double sumR = 0;
    for(int32_t j = 0; j < nComp; j++) {
      double sigmaSq = SIGMASQ(comps[j]);
      double a = y - MU(comps[j]);
      a = 4*sigmaSq + a*a;
      a = 1./a;
      double b = sigmaSq*a;
      double Rj = comps[j].w*b*b*sqrt(a);
      sumR += Rj;
      R[j] = Rj;
      Ruy[j] = 5*a*Rj;
    }
    sumR = 1./sumR;
    for(int32_t j = 0; j < nComp; j++) {
      chi[j] += R[j]*sumR;
      double Ruyj = Ruy[j]*sumR;
      n[j] += Ruyj;
      s[j] += y * Ruyj;
    }
  }
  *phi = nY;

  if(U->s == U->e)
    return;

  double P0 = 1;
  double H3[nComp], H1[nComp], H0[nComp];
  memset(H3, 0, sizeof H3);
  memset(H1, 0, sizeof H1);
  memset(H0, 0, sizeof H0);
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

      H0[j] += cst*(fun0(bSin) - fun0(aSin));
      H1[j] += (fun1(aNorm,4) - fun1(bNorm,4)) * cst / 5;
    }
  }
  *phi = nY / P0;
  
  for(int32_t j = 0; j < nComp; j++) {
    double phiW = (*phi)*comps[j].w;
    chi[j] += phiW * H3[j]; 
    s[j] += phiW * (MU(comps[j])*H0[j] + 2*sqrt(SIGMASQ(comps[j]))*H1[j]) / SIGMASQ(comps[j]);
    n[j] += phiW * H0[j] / SIGMASQ(comps[j]);
  }
}
