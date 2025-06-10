#include <cmath>
#include <cstdint>
#include <omp.h>
#include <valarray>

extern float cnv;

void SFR_Buffered(std::valarray<float> &sfr, uint32_t buff_sz, const float z) {

  uint_fast64_t i;

  const float a = 2.59, b = -0.62, d = 0.4, e = -2.25, c1 = 8e8, c2 = 7e9,
              c3 = 1e11;

  const float fac = 2.25e-26 * (1.0 + 7.5e-2 * (z - 7.0));

#pragma omp parallel num_threads(4)
  {

#pragma omp for
    for (i = 0; i < buff_sz; ++i) {

      sfr[i] *= cnv;
      sfr[i] = powf(sfr[i], a) * powf(1 + sfr[i] / c1, b) *
               powf(1 + sfr[i] / c2, d) * powf(1 + sfr[i] / c3, e) * fac;
    }
  }

} // End of SFR_Buffered()
