#include "infmat_miss.c"

static int getInfMat(struct InputData* seg, struct MixtureComp* comps, int32_t nComp, double* infMat)
{
  int32_t NF = seg->P.e - seg->P.s;
  int32_t NR = seg->N.e - seg->N.s;

  int32_t dim = 5*nComp-1;
  memset(infMat, 0, sizeof(double)*dim*dim);

  double scoreObsBarF[dim];
  memset(scoreObsBarF, 0, sizeof(double)*dim);

  int r,c;
  int i,j;
  for(i = 0; i < NF; i++) {
    double rF[nComp];
    double ruyF[nComp];
    double qWoF[nComp];
    double qMuoF[nComp];
    double qSigmaoF[nComp];
    double rFSum = 0;
    for(j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double fMu = c->mu - c->delta/2;
      double sigma = sqrt(c->sigmaSqF);
      double yNormF = (seg->P.s[i] - fMu)/sigma;
      rF[j] = c->w * tdist4(yNormF) /sigma;
      rFSum += rF[j];
    }
    for(j = 0; j < nComp; j++) {
      rF[j] /= rFSum;
      struct MixtureComp* c = comps + j;
      double fMu = c->mu - c->delta/2;
      double sigma = sqrt(c->sigmaSqF);
      double yNormF = (seg->P.s[i] - fMu)/sigma;
      ruyF[j] = rF[j] * 5.0/(4 + yNormF*yNormF);
      qWoF[j] = rF[j] / c->w;
      qMuoF[j] = ruyF[j]*yNormF/sigma;
      qSigmaoF[j] = c->sigmaSqF/2.*(rF[j]-ruyF[j]*yNormF*yNormF);
    }

    double scoreObsF[dim];
    memset(scoreObsF, 0, sizeof scoreObsF);
    for(j=1;j<nComp;j++)
      scoreObsF[j-1] = qWoF[j];
    for(j=0;j<nComp;j++)
    {
      scoreObsF[nComp-1+j] = qMuoF[j];
      scoreObsF[2*nComp-1+j] = -qMuoF[j]/2.;
      scoreObsF[3*nComp-1+j] = qSigmaoF[j];
    }

    for(r=0;r<dim;r++)
       scoreObsBarF[r] += scoreObsF[r];
    for(r=0;r<dim;r++)
      for(c=0;c<=r;c++)
        infMat[r*dim+c] += scoreObsF[r]*scoreObsF[c];
  }

  for(r=0;r<dim;r++)
    for(c=0;c<=r;c++)
      infMat[r*dim+c] -= scoreObsBarF[r]*scoreObsBarF[c]/NF;

  double scoreObsBarR[dim];
  memset(scoreObsBarR, 0, sizeof(double)*dim);

  for(i = 0; i < NR; i++) {
    double rR[nComp];
    double ruyR[nComp];
    double qWoR[nComp];
    double qMuoR[nComp];
    double qSigmaoR[nComp];
    double rRSum = 0;
    for(j = 0; j < nComp; j++) {
      struct MixtureComp* c = comps + j;
      double rMu = c->mu + c->delta/2;
      double sigma = sqrt(c->sigmaSqR);
      double yNormR = (seg->N.s[i] - rMu)/sigma;
      rR[j] = c->w * tdist4(yNormR) /sigma;
      rRSum += rR[j];
    }
    for(j = 0; j < nComp; j++) {
      rR[j] /= rRSum;
      struct MixtureComp* c = comps + j;
      double rMu = c->mu + c->delta/2;
      double sigma = sqrt(c->sigmaSqR);
      double yNormR = (seg->N.s[i] - rMu)/sigma;
      ruyR[j] = rR[j] * 5.0/(4 + yNormR*yNormR);
      qWoR[j] = rR[j] / c->w;
      qMuoR[j] = ruyR[j]*yNormR/sigma;
      qSigmaoR[j] = c->sigmaSqR/2.*(rR[j]-ruyR[j]*yNormR*yNormR);
    }

    double scoreObsR[dim];
    memset(scoreObsR, 0, sizeof scoreObsR);
    for(j=1;j<nComp;j++)
      scoreObsR[j-1] = qWoR[j];
    for(j=0;j<nComp;j++)
    {
      scoreObsR[nComp-1+j] = qMuoR[j];
      scoreObsR[2*nComp-1+j] = qMuoR[j]/2.;
      scoreObsR[4*nComp-1+j] = qSigmaoR[j];
    }

    for(r=0;r<dim;r++)
       scoreObsBarR[r] += scoreObsR[r];
    for(r=0;r<dim;r++)
      for(c=0;c<=r;c++)
        infMat[r*dim+c] += scoreObsR[r]*scoreObsR[c];
  }

  for(r=0;r<dim;r++)
    for(c=0;c<=r;c++)
      infMat[r*dim+c] -= scoreObsBarR[r]*scoreObsBarR[c]/NR;
 
  double infMatPri[dim*dim];
  memset(infMatPri, 0, sizeof(double)*dim*dim);
  for(j=0;j<nComp;j++)
  {
    struct MixtureComp c = comps[j];
    infMatPri[(2*nComp-1+j)*dim+2*nComp-1+j] = PriorParams.rho*(1./c.sigmaSqF+1./c.sigmaSqR);
    infMatPri[(3*nComp-1+j)*dim+3*nComp-1+j] = (PriorParams.alpha-1./2.)*c.sigmaSqF*c.sigmaSqF;
    infMatPri[(4*nComp-1+j)*dim+4*nComp-1+j] = (PriorParams.alpha-1./2.)*c.sigmaSqR*c.sigmaSqR;
    infMatPri[(2*nComp-1+j)*dim+3*nComp-1+j] = (c.delta-PriorParams.xi)*PriorParams.rho;
    infMatPri[(2*nComp-1+j)*dim+4*nComp-1+j] = (c.delta-PriorParams.xi)*PriorParams.rho;
    infMatPri[(3*nComp-1+j)*dim+2*nComp-1+j] = (c.delta-PriorParams.xi)*PriorParams.rho;
    infMatPri[(4*nComp-1+j)*dim+2*nComp-1+j] = (c.delta-PriorParams.xi)*PriorParams.rho;
  }
  if(nComp>1)
  {
    infMatPri[(nComp)*dim+nComp] = 0;
    infMatPri[(2*nComp-1)*dim+2*nComp-1] = 0;
    
    for(j=1;j<(nComp-1);j++)
      infMatPri[(nComp+j)*dim+nComp+j] = 0;

    for(j=0;j<(nComp-1);j++)
    {
      infMatPri[(nComp+j+1)*dim+nComp+j] = 0;
      infMatPri[(nComp+j)*dim+nComp+j+1] = 0;
    }
  }
  for(r=0;r<dim;r++)
    for(c=0;c<dim;c++)
      infMat[r*dim+c] += infMatPri[r*dim+c];

  if(seg->U.s != seg->U.e) {
    addInfMatMissF(seg, comps, nComp, infMat);
    addInfMatMissR(seg, comps, nComp, infMat);
  }
  
  gsl_matrix_view gslInfMat = gsl_matrix_view_array(infMat, dim, dim);
  int flag = gsl_linalg_cholesky_decomp(&gslInfMat.matrix);
  if(flag==GSL_EDOM)
    return(flag);

  gsl_matrix *M=gsl_matrix_calloc(dim, dim);
  gsl_matrix_set_identity(M);
  /* Compute the inverse of I^{-1/2} */
  /* Compute L'^{-1} */
  flag=gsl_blas_dtrsm(CblasLeft, CblasUpper, CblasNoTrans, CblasNonUnit, 1.0, &gslInfMat.matrix, M);
  if(flag!=0)
    return flag;

  gsl_vector *A=gsl_vector_calloc(dim);
  for(j=0;j<nComp;j++)
  {
    struct MixtureComp* c = comps+j;
    gsl_vector_set_zero(A);
    gsl_vector_set(A,nComp-1+j,1);
    /* Compute L'^{-1} A */
    flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, M, A);
    if(flag!=0)
      return flag;
    /* Compute the norm of that vector */
    c->se = gsl_blas_dnrm2(A);
    // Make sure it's a number
    if(gsl_finite(c->se)==0)
      return 1;

    gsl_vector_set_zero(A);
    gsl_vector_set(A,nComp-1+j,1);
    gsl_vector_set(A,nComp-1+nComp+j,-0.5);
    flag=gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, M, A);
    if(flag!=0)
      return flag;
    c->seF = gsl_blas_dnrm2(A);
    if(gsl_finite(c->seF)==0)
      return 1;

    gsl_vector_set_zero(A);
    gsl_vector_set(A,nComp-1+j,1);
    gsl_vector_set(A,nComp-1+nComp+j,0.5);
    flag=gsl_blas_dtrmv (CblasUpper, CblasNoTrans, CblasNonUnit, M, A);
    if(flag!=0)
      return flag;
    c->seR = gsl_blas_dnrm2(A);
    if(gsl_finite(c->seR)==0)
      return 1;
  }
  
  // Copy infMat in M
  // So that I don't have to recompute the inverse of L'
  gsl_matrix_memcpy(&gslInfMat.matrix, M);
  
  gsl_vector_free(A);
  gsl_matrix_free(M);

  return(flag);
}
