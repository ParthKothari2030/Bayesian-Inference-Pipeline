#include "pwr_spc_decl.H"

moments::moments () {}

moments::moments (unsigned N)
{
  k_mode.resize (N);
  k_mode = 0;
  P0.resize (N);
  P0 = 0;
  P2.resize (N);
  P2 = 0;
  P4.resize (N);
  P4 = 0;
  no.resize (N);
  no = 0;
}

void
moments::write (std::string filename)
{
  std::ofstream out;
  out.open (filename);

  out.fill ('_');
  out << std::left << std::setw (12) << "#k (1/Mpc)"
      << "  " << std::setw (12) << "P0 "
      << "  " << std::setw (12) << "P2 "
      << "  " << std::setw (12) << "P4 "
      << "  " << std::setw (12) << "N " << '\n';
  out.fill (' ');

  for (unsigned i = 0; i < k_mode.size (); ++i)
    {
      out << std::uppercase << std::left << std::showpos << std::setw (12)
          << k_mode[i] << "  " << std::setw (12) << P0[i] << "  "
          << std::setw (12) << P2[i] << "  " << std::setw (12) << P4[i]
          << "  ";
      out << std::setw (12) << no[i] << "\n";
    }

  out.close ();
}

void
moments::clear ()
{
  k_mode.resize (0);
  P0.resize (0);
  P2.resize (0);
  P4.resize (0);
  no.resize (0);
}

moments::~moments () {}
