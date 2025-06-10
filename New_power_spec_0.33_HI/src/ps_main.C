#include "ps_aux.H"
#include "pwr_spc.H"

int N1, N2, N3;
float l, k_min, k_max;
unsigned Nbin;

int
main (int argc, char *argv[])
{
  std::string dir, command, path, filename, outdir;

  std::ifstream in{ argv[1], std::ios_base::in };

  if (in.is_open ()==false)
    {
      std::cout << "Error! Input file not given"
                   "...terminating program.\n";
      exit (1);
    }

  in >> Nbin >> k_min >> k_max;

  val_f box;

  int N[3];
  int counter = 0;

  moments pw_mom;

  while (in >> path >> filename >> outdir)
    {
      box = read_box (path + filename);
      subtract_mean (box);
      N[0] = N1;
      N[1] = N2;
      N[2] = N3;

      // std::cout << "k_min = " << k_min << "\nk_max = " << k_max;
      // std::cout << "\nN_bin = " << Nbin << std::endl;
      // std::cout << N1 << "\t" << N2 << "\t" << N3 << "\t" << l << std::endl;

      pw_mom = cal_pow_updated (box, N, l, k_min, k_max);
      command = "mkdir -p " + outdir;
      size_t bytt = system (command.c_str ());
      pw_mom.write (outdir + '/' + filename + ".powspec");
      pw_mom.clear ();

      ++counter;
      // std::cout << "No. of files done: " << counter << '\n';
    } // end while

  in.close ();

  // std::cout << "chain Power spectrum done\n";
} // End of main ()
