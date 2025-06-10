	#ifndef PS_AUX_H_
	#define PS_AUX_H_

	#include<iostream>
	#include<fstream>
	#include<string>
	#include<valarray>
	#include<fftw3.h>

	using val_f =std::valarray<float>;

	val_f read_box(const std::string);
	void subtract_mean(val_f&);

	#endif
