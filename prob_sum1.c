static inline void FUNCTION_NAME (
  struct MixtureComp* comps, int32_t nComp,
  int32_t* Y, int32_t nY,
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
}
