static void addInfMatMissF(struct InputData* seg, struct MixtureComp* comps, int32_t nComp,
 double* infMat) {
  int32_t J = seg->U.e - seg->U.s;
  int32_t dim = 5*nComp-1;
  int32_t nY = seg->P.e - seg->P.s;

  double infMatMiss[dim*dim];
  memset(infMatMiss, 0, sizeof(double)*dim*dim);

  double scoreMissBar[dim];
  memset(scoreMissBar, 0, sizeof(double)*dim);

  double P0 = 1;
  for(struct Interval* u = seg->U.s; u < seg->U.e; u++) {
    int32_t a = max(u->start, seg->bounds.start);
    int32_t b = min(u->end, seg->bounds.end) - 1;
    double PJ = 0;
    double H0[nComp], H1[nComp], H2[nComp], H3[nComp];
    memset(H0, 0, sizeof H0);
    memset(H1, 0, sizeof H1);
    memset(H2, 0, sizeof H2);
    memset(H3, 0, sizeof H3);
    for(int32_t j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double sigma = sqrt(c->sigmaSqF);
      double mu = c->mu - 0.5*c->delta;
      double aNorm = (a - mu) / sigma;
      double bNorm = (b - mu) / sigma;
      double aSin = sin(atan(aNorm/2));
      double bSin = sin(atan(bNorm/2));
      H0[j] = cst*(fun0(bSin) - fun0(aSin));
      H1[j] = (fun1(aNorm,4) - fun1(bNorm,4)) * cst / 5;
      H2[j] = cst*(fun2(bSin) - fun2(aSin));
      H3[j] = gsl_cdf_tdist_P(bNorm, 4) - gsl_cdf_tdist_P(aNorm, 4);
      PJ += c->w * H3[j];
    }
    P0 -= PJ;

    double scoreMiss[dim];
    for(int32_t j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double sigma = sqrt(c->sigmaSqF);
      double wP = c->w / PJ;
      double qWm = H3[j] / PJ;
      double qMum = wP * H1[j] * 2. / sigma;
      double qSigmam = wP * c->sigmaSqF/2.*(H0[j] - 4.*H2[j]);
      if(j > 0)
        scoreMiss[j-1] = qWm;
      scoreMiss[nComp-1+j] = qMum;
      scoreMiss[2*nComp-1+j] = - qMum/2.;
      scoreMiss[3*nComp-1+j] = qSigmam;
    }

    for(int32_t i = 0; i < dim; i++)
      for(int32_t j = 0; j <= i; j++)
        infMatMiss[i*dim+j] += PJ*nY*scoreMiss[i]*scoreMiss[j];

    for(int32_t i = 0; i < dim; i++)
      scoreMissBar[i] += PJ*scoreMiss[i]/J;
  }

  if(P0 == 1)
    return;

  for(int32_t i = 0; i < dim; i++)
    for(int32_t j = 0; j <= i; j++)
      infMatMiss[i*dim+j] /= P0;

  for(int32_t i = 0; i < dim; i++)
    for(int32_t j = 0; j <= i; j++)
      infMatMiss[i*dim+j] -= nY*scoreMissBar[i]*scoreMissBar[j]/P0/(1-P0);
    
  for(int32_t i = 0; i < dim; i++)
    for(int32_t j = 0; j <= i; j++)
      infMat[i*dim+j] += infMatMiss[i*dim+j];
}


static void addInfMatMissR(struct InputData* seg, struct MixtureComp* comps, int32_t nComp,
 double* infMat) {
  int32_t J = seg->U.e - seg->U.s;
  int32_t dim = 5*nComp-1;
  int32_t nY = seg->N.e - seg->N.s;

  double infMatMiss[dim*dim];
  memset(infMatMiss, 0, sizeof(double)*dim*dim);

  double scoreMissBar[dim];
  memset(scoreMissBar, 0, sizeof(double)*dim);

  double P0 = 1;
  for(struct Interval* u = seg->U.s; u < seg->U.e; u++) {
    int32_t a = max(u->start, seg->bounds.start);
    int32_t b = min(u->end, seg->bounds.end) - 1;
    double PJ = 0;
    double H0[nComp], H1[nComp], H2[nComp], H3[nComp];
    memset(H0, 0, sizeof H0);
    memset(H1, 0, sizeof H1);
    memset(H2, 0, sizeof H2);
    memset(H3, 0, sizeof H3);
    for(int32_t j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double sigma = sqrt(c->sigmaSqR);
      double mu = c->mu + 0.5*c->delta;
      double aNorm = (a - mu) / sigma;
      double bNorm = (b - mu) / sigma;
      double aSin = sin(atan(aNorm/2));
      double bSin = sin(atan(bNorm/2));
      H0[j] = cst*(fun0(bSin) - fun0(aSin));
      H1[j] = (fun1(aNorm,4) - fun1(bNorm,4)) * cst / 5;
      H2[j] = cst*(fun2(bSin) - fun2(aSin));
      H3[j] = gsl_cdf_tdist_P(bNorm, 4) - gsl_cdf_tdist_P(aNorm, 4);
      PJ += c->w * H3[j];
    }
    P0 -= PJ;

    double scoreMiss[dim];
    for(int32_t j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double sigma = sqrt(c->sigmaSqR);
      double wP = c->w / PJ;
      double qWm = H3[j] / PJ;
      double qMum = wP * H1[j] * 2. / sigma;
      double qSigmam = wP * c->sigmaSqR/2.*(H0[j] - 4.*H2[j]);
      if(j > 0)
        scoreMiss[j-1] = qWm;
      scoreMiss[nComp-1+j] = qMum;
      scoreMiss[2*nComp-1+j] = qMum/2.;
      scoreMiss[4*nComp-1+j] = qSigmam;
    }

    for(int32_t i = 0; i < dim; i++)
      for(int32_t j = 0; j <= i; j++)
        infMatMiss[i*dim+j] += PJ*nY*scoreMiss[i]*scoreMiss[j];

    for(int32_t i = 0; i < dim; i++)
      scoreMissBar[i] += PJ*scoreMiss[i]/J;
  }

  if(P0 == 1)
    return;

  for(int32_t i = 0; i < dim; i++)
    for(int32_t j = 0; j <= i; j++)
      infMatMiss[i*dim+j] /= P0;

  for(int32_t i = 0; i < dim; i++)
    for(int32_t j = 0; j <= i; j++)
      infMatMiss[i*dim+j] -= nY*scoreMissBar[i]*scoreMissBar[j]/P0/(1-P0);
    
  for(int32_t i = 0; i < dim; i++)
    for(int32_t j = 0; j <= i; j++)
      infMat[i*dim+j] += infMatMiss[i*dim+j];
}
