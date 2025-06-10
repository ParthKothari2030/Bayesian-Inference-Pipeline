	#ifndef PWR_SPC_DECL_H_
	#define PWR_SPC_DECL_H_

	#include<valarray>
	#include<string>
	#include<fstream>
	#include<iomanip>

	struct moments {
	   std::valarray<double>
	   k_mode, P0, P2, P4;
	   std::valarray<std::size_t> no;

	   moments();
	   moments(unsigned N);

	   void write(std::string filename);
	   void clear();
	   ~moments();
	};

	#endif
