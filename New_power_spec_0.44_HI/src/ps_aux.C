#include <cstdint>
#include <numeric>
#include "ps_aux.H"

extern int N1, N2, N3;
extern float l;

inline long double
add (long double x, float& y)
{
  return x + y;
}

inline uint64_t
fbt (uint64_t Num)
{
  return Num * sizeof (float);
}

val_f
read_box (const std::string file)
{
  std::ifstream in{ file, std::ios_base::binary };

  if (!in.is_open ())
    {
      std::cout << "Error: File " << file
                << " not found...terminating program.\n";
      exit (1);
    }

  uint64_t N[3];
  in.read ((char *)N, sizeof (uint64_t) * 3);

  N1 = N[0];
  N2 = N[1];
  N3 = N[2];

  val_f box (N1 * N2 * N3);

  in.read ((char *)&l, sizeof (float));

  in.read ((char *)&box[0], sizeof (float) * N[0] * N[1] * N[2]);
  in.close ();

  return box;
}

void
subtract_mean (val_f &box)
{
  uint64_t num;
  uint64_t i;
  uint64_t N[3];

  N[0] = N1;
  N[1] = N2;
  N[2] = N3;

  long double sum_mean;

  std::valarray<float> slice;

  for (num = 0; num < N3; ++num)
    {
      slice = box[std::gslice (num, { N[0] * N[1] }, { N[2] })];

      sum_mean = std::accumulate (begin (slice), end (slice), 0., add) / slice.size ();
//      for (i = 0, sum_mean = 0; i < slice.size (); ++i)
//        sum_mean += slice[i];

//      sum_mean /= slice.size ();
      slice -= sum_mean;

      box[std::gslice (num, { N[0] * N[1] }, { N[2] })] = slice;
    } // end for
} // end of subtract_mean()
