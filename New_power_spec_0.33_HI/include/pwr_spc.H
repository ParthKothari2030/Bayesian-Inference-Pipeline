	#ifndef PWR_SPC_H_
	#define PWR_SPC_H_

	#include<iostream>
	#include<cmath>
	#include<valarray>
	#include<complex>
	#include<algorithm>
	#include<fftw3.h>
	#include<omp.h>
	#include"pwr_spc_decl.H"

	using val_f =std::valarray<float>;

	moments cal_pow_updated(val_f&, const int [], const float, float, float);

	#endif
