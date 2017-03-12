static inline void FUNCTION_NAME (
  struct MixtureComp* comps, int32_t nComp,
  int32_t* Y, int32_t nY,
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
}
