int mergePeaks(struct MixtureComp* comps, int32_t* nComp, double* infMat)
{
  const int nSe = 3;

  int32_t K = *nComp;
  int32_t K0 = K;

  double I[K*K];
  memset(I, 0, sizeof(double)*K*K);
  int i,j;
  for(i=0;i<K;i++)
    I[i*K+i] = 1;

  int dim = 5*K - 1;

  gsl_matrix_view gslInfMat = gsl_matrix_view_array(infMat, dim, dim);

  double A[dim], B[dim], C[dim];
  gsl_vector_view gslA = gsl_vector_view_array(A, dim);
  gsl_vector_view gslB = gsl_vector_view_array(B, dim);
  gsl_vector_view gslC = gsl_vector_view_array(C, dim);

  double origW[K];
  for(i = 0;i < K; i++)
    origW[i] = comps[i].w;

  for(;;) {
    double minDiff = 0;
    int iMerge = -1;
    for(i = 0; i < K - 1; i++) {
      struct MixtureComp cur = comps[i];
      struct MixtureComp next = comps[i+1];
      double diff = (next.mu - next.delta/2. - nSe*next.seF) - (cur.mu + cur.delta/2. + nSe*cur.seR);
      if(diff < minDiff) {
        minDiff = diff;
        iMerge = i;
      }
    }

    if(minDiff >= 0)
      break;

    // Compute the new parameters
    struct MixtureComp* cur = comps + iMerge;
    struct MixtureComp* next =  cur + 1;
    double SumCombW = cur->w + next->w;
    double TmpMu = (cur->w*cur->mu + next->w*next->mu)/SumCombW;
    double TmpDelta=(cur->w*cur->delta + next->w*next->delta)/SumCombW;
    double nu = 4;
    double TmpSigmaSqF=(nu-2.)/nu*(((nu/(nu-2.)*cur->sigmaSqF+gsl_pow_2(cur->mu-cur->delta/2.))*cur->w+(nu/(nu-2.)*next->sigmaSqF+gsl_pow_2(next->mu-next->delta/2.))*next->w)/SumCombW-gsl_pow_2(TmpMu-TmpDelta/2.));
   // double TmpSigmaSqF=0.5*(((2.*cur->sigmaSqF+gsl_pow_2(cur->mu-cur->delta/2.))*cur->w+(2.*next->sigmaSqF+gsl_pow_2(next->mu-next->delta/2.))*next->w)/SumCombW-gsl_pow_2(TmpMu-TmpDelta/2.));
    double TmpSigmaSqR=0.5*(((2.*cur->sigmaSqR+gsl_pow_2(cur->mu+cur->delta/2.))*cur->w+(2.*next->sigmaSqR+gsl_pow_2(next->mu+next->delta/2.))*next->w)/SumCombW-gsl_pow_2(TmpMu+TmpDelta/2.));
    cur->mu=TmpMu;
    cur->delta=TmpDelta;
    cur->w=SumCombW;
    cur->sigmaSqF=TmpSigmaSqF;
    cur->sigmaSqR=TmpSigmaSqR;
    
    double SumW=0;
    for(i=0;i<K0;i++)
      if(I[iMerge*K0+i] || I[(iMerge+1)*K0+i])
      {
        I[iMerge*K0+i] = 1;
        SumW += origW[i];
      }

    memset(A, 0, sizeof(double)*dim);
    memset(B, 0, sizeof(double)*dim);
    memset(C, 0, sizeof(double)*dim);
    for(i=0;i<K0;i++)
      if(I[iMerge*dim+i])
      {
        A[K0-1+i] = origW[i]/SumW;
        B[K0-1+i] = origW[i]/SumW;
        B[K0-1+K0+i] = -origW[i]/SumW/2.;
        C[K0-1+i] = origW[i]/SumW;
        C[K0-1+K0+i] = origW[i]/SumW/2.;
      }

   //printf("A:\n");
   //for(i=0;i<dim;i++)
   //  printf(" %.4f", A[i]);
   //printf("\n");

   //printf("infMat:\n");
   //for(i=0;i<dim;i++) {
   //  for(j=0;j<dim;j++)
   //    printf(" %.4f", infMat[i*dim+j]);
   //  printf("\n");
   //}

    int flag = gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &gslInfMat.matrix, &gslA.vector);
    if(flag!=0)
      return flag;
    cur->se = gsl_blas_dnrm2(&gslA.vector);

    flag = gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &gslInfMat.matrix, &gslB.vector);
    if(flag!=0)
      return flag;
    cur->seF = gsl_blas_dnrm2(&gslB.vector);

    flag = gsl_blas_dtrmv(CblasUpper, CblasNoTrans, CblasNonUnit, &gslInfMat.matrix, &gslC.vector);
    if(flag!=0)
      return flag;
    cur->seR = gsl_blas_dnrm2(&gslC.vector);

    /* Shift the components */
    for(i=(iMerge+1);i<K-1;i++)
    {
      comps[i] = comps[i+1];
      for(j=0;j<K;j++)// !!! j<K0
        I[i*K0+j] = I[(i+1)*K0+j];
    }

    /* Decrease the number of components by one */
    K--;               
  }

  *nComp = K;
  
  return 0;
}
