#ifndef _READ_HALO_H
#define _READ_HALO_H

#include <string>
#include <valarray>
#include    <array>

void File_Handler( std::array<std::ifstream,4>&,
									 int8_t                       ) ;

void  Read_Halos_Buffered( std::array<std::ifstream,4>&,
													 std::array<char*,4>&,
													 uint32_t                     ) ;
#endif
