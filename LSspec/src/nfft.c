#include <math.h>
#include <complex.h>
#include <nfft3.h>
 

void nfft(const double* t, const double* y, int* N, int* M, complex double* d)
{
  /* t are the times
   * y are the time-domain values
   * N is the number of Fourier coefs
   * M is the number of non-equispaced nodes
   * d are the frequency-domain values
   */
  nfft_plan p;
  nfft_init_1d(&p, *N, *M);

  // set the node locations (equispaced)
  for (int i = 0; i < *M; i++)
      p.x[i] = t[i];
  // set the input values
  for (int i = 0; i < *M; i++)
      p.f[i] = y[i];

  //if (p.nfft_flags & PRE_ONE_PSI)
    nfft_precompute_one_psi(&p);
  
  // compute the nfft
  nfft_adjoint(&p);
  
  // retrieve the nfft
  for (int i = 0; i < *N; i++)
    d[i] = p.f_hat[i];

  nfft_finalize(&p);
}

void nfft2(const double* t, const double* y, int* N, int* M, complex double* d)
{
  /* t are the times
  * y are the time-domain values
  * N is the number of Fourier coefs
  * M is the number of non-equispaced nodes
  * d are the frequency-domain values
  */
  nfft_plan p;
  nfft_init_1d(&p, *N, *M);
  
  // set the node locations (equispaced)
  for (int i = 0; i < *M; i++)
    p.x[i] = t[i];
  // set the input values
  for (int i = 0; i < *M; i++)
    p.f_hat[i] = y[i];
  
  //if (p.nfft_flags & PRE_ONE_PSI)
  nfft_precompute_one_psi(&p);
  
  // compute the nfft
  nfft_trafo(&p);
  
  // retrieve the nfft
  for (int i = 0; i < *N; i++)
    d[i] = p.f[i];
  
  nfft_finalize(&p);
}
