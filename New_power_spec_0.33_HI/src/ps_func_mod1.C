#include "pwr_spc.H"

extern unsigned Nbin;

inline float
sqr (float x)
{
  return x * x;
};

inline float
cube (float x)
{
  return x * x * x;
};

moments
cal_pow_updated (val_f &ro_dum, const int N[], const float LL, float k_min, float k_max)
{

  std::valarray<std::complex<float> > ro_dum_fourier (
      N[0] * N[1] * (N[2] / 2 + 1)); // allocate memory for output

  int q = fftwf_init_threads ();
  if (q != 0)
    fftwf_plan_with_nthreads (omp_get_max_threads ());

  // TAKING FOURIER TRANSFORM OF ro_dum
  fftwf_plan p_ro_dum = fftwf_plan_dft_r2c (
      3, N, &ro_dum[0], (fftwf_complex *)&ro_dum_fourier[0], FFTW_ESTIMATE);

  fftwf_execute (p_ro_dum);
  ro_dum.resize (0);
  fftwf_destroy_plan (p_ro_dum);

  moments pw_mom (Nbin);

  std::valarray<double> no2 (0., Nbin);
  std::valarray<double> no4 (0., Nbin);

  int N_max = std::max ({ N[0], N[1], N[2] });
  
  float R_log, logspace_bin_w;
  
  k_min *= LL / (2 * M_PI);
  k_max *= LL / (2 * M_PI);
  // std::cout<<"LL= "<<LL<<"\t"<<"\n"<<"k_min_grid= "<<k_min<<"\n k_max_grid= "<<k_max;
  
  if (k_min <= 1./N_max || k_max >= 1./2){
    std::cout<<"Invalid k values.\n";
    exit(1);
  }   

  R_log = log10 (k_max) - log10 (k_min);
  logspace_bin_w = R_log / Nbin;

  // //---------------------- BINNING POWER SPECTRA -------------------//

  float kk, mu, P0, P2, P4;
  long index, i, j, k, a, b, c, bin;

  // half-lines
  for (i = 1; i <= N[0] / 2; ++i)
    for (j = 0; j <= N[1] / 2; j += N[1] / 2)
      for (k = 0; k <= N[2] / 2; k += N[2] / 2)
        {

          kk = sqrt (sqr (1. * i / N[0]) + sqr (1. * j / N[1])
                     + sqr (1. * k / N[2]));
          if (kk < k_min || kk > k_max)
            continue;

          index = k + (N[2] / 2 + 1) * (j + N[1] * i);
          
          mu = (1. * k) / (N[2] * kk);

          P2 = 0.5 * (3. * mu * mu - 1.);
          P4 = 0.125 * (35. * powf (mu, 4.) - 30. + 0 * mu * mu + 3);

          bin = (log10 (kk) - log10 (k_min)) / logspace_bin_w;
          // logarithmic bins
          // std::cout<<"Bin= "<<bin<<"\n";

          P0 = (double)norm (ro_dum_fourier[index]);
          if(bin <= Nbin - 1 )
          {
            pw_mom.P0[bin] += P0;
            pw_mom.P2[bin] += P0 * (double)P2;
            pw_mom.P4[bin] += P0 * (double)P4;

            pw_mom.k_mode[bin] += (double)kk;
            // std::cout <<"kmode_halfline= "<< pw_mom.k_mode[bin] << "\n";
            pw_mom.no[bin] += 1;

            no2[bin] += (double)(P2 * P2);
            no4[bin] += (double)(P4 * P4);
          }
        } // end of half-lines

  // half-planes
  for (i = 0; i < N[0]; ++i)
    {

      a = (i > N[0] / 2) ? N[0] - i : i;
      for (j = 1; j < N[1] / 2; ++j)
        {

          b = j;
          for (k = 0; k <= N[2] / 2; k += N[2] / 2)
            {

              c = k;
              kk = sqrt (sqr (1. * a / N[0]) + sqr (1. * b / N[1])
                         + sqr (1. * c / N[2]));
              if (kk < k_min || kk > k_max)
                continue;

              index = k + (N[2] / 2 + 1) * (j + N[1] * i);
              
              mu = (1. * c) / (N[2] * kk);

              P2 = 0.5 * (3. * mu * mu - 1.);
              P4 = 0.125 * (35. * powf (mu, 4.) - 30. * mu * mu + 3);

              bin = (log10 (kk) - log10 (k_min)) / logspace_bin_w;
              // logarithmic bins
              // std::cout<<"\n K= "<<kk;
              P0 = (double)norm (ro_dum_fourier[index]);
              if(bin <= Nbin-1)
              {
                pw_mom.P0[bin] += P0;
                pw_mom.P2[bin] += P0 * (double)P2;
                pw_mom.P4[bin] += P0 * (double)P4;

                pw_mom.k_mode[bin] += (double)kk;
                // std::cout <<"kmode_halfplane= "<< pw_mom.k_mode[bin] << "\n";
                pw_mom.no[bin] += 1;
                // std::cout<<"\n number= "<<pw_mom.no[bin]<<"\t";

                no2[bin] += (double)(P2 * P2);
                no4[bin] += (double)(P4 * P4);
              }
            }
        }
    } // end of half-planes

  // half-cube
  for (i = 0; i < N[0]; ++i)
    {

      a = (i > N[0] / 2) ? N[0] - i : i;
      for (j = 0; j < N[1]; ++j)
        {

          b = (j > N[1] / 2) ? N[1] - j : j;
          for (k = 1; k < N[2] / 2; ++k)
            {

              c = k;
              kk = sqrt (sqr (1. * a / N[0]) + sqr (1. * b / N[1])
                         + sqr (1. * c / N[2]));
              if (kk < k_min || kk > k_max)
                continue;

              index = k + (N[2] / 2 + 1) * (j + N[1] * i);
              // std::cout<<"\n kk= "<<kk;
              // kk*(2 * pi/LL) is |k|
              // kk=1/2 corresponds to kmode[Nbin-1] i.e. Nyquist

              mu = (1.0 * c) / (N[2] * kk);

              P2 = 0.5 * (3. * mu * mu - 1.);
              P4 = 0.125 * (35.0 * powf (mu, 4.0) - 30.0 * mu * mu + 3);

              bin = (log10 (kk) - log10 (k_min)) / logspace_bin_w;
              // logarithmic bins

              P0 = (double)norm (ro_dum_fourier[index]);
              if(bin <= Nbin - 1)
              {
                pw_mom.P0[bin] += P0;
                // std::cout<<"\n P0= "<<pw_mom.P0[bin]<<"\t";
                // std::cout<<"\n bin = "<<bin<<"\t";
                pw_mom.P2[bin] += P0 * (double)P2;
                pw_mom.P4[bin] += P0 * (double)P4;

                pw_mom.k_mode[bin] += (double)kk;
                // std::cout <<"kmode_halfcube= "<< pw_mom.k_mode[bin] << "\n";
                pw_mom.no[bin] += 1;
                // std::cout<<"\n number= "<<pw_mom.no[bin]<<"\t";

                no2[bin] += (double)(P2 * P2);
                no4[bin] += (double)(P4 * P4);
              }
            }
        }
    } // end of half-cubes

  double norml = 1. / (N[0] * N[1] * N[2]);
  for (i = 0; i < Nbin; ++i)
    {
      if (pw_mom.no[i] > 0)
        {
          pw_mom.P0[i] *= cube (LL) * norml / pw_mom.no[i];
          // std::cout <<"pw_mom ="<< pw_mom.no[i] << "\n";
          // std::cout<<"\n k_mode ="<<pw_mom.k_mode[i];
          
          pw_mom.k_mode[i] *= 2 * M_PI / (pw_mom.no[i] * LL);
          // std::cout<<"k_mode_MPc ="<<pw_mom.k_mode[i]<<"\n";
        }

      if (no2[i] > 0.0)
        pw_mom.P2[i] *= cube (LL) * norml / (no2[i]);

      if (no4[i] > 0.0)
        pw_mom.P4[i] *= cube (LL) * norml / (no4[i]);
    }

  return pw_mom;
} // end of cal_pow_updated()