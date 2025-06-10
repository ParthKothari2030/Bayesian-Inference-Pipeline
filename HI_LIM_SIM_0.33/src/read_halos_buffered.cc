#include <array>
#include <cstdint>
#include <fstream>
#include <omp.h>
#include <string>

extern std::string Dir;

char Buffer[4][4096];

void File_Handler (std::array<std::ifstream, 4> &in, int8_t flag) {

  std::string File[] = {"Pos_x.bin", "Pos_y.bin", "Pos_z.bin", "Mass.bin"};

  if (flag == 0) {

    for (uint8_t i = 0; i < 4; ++i) {
      in[i].rdbuf()->pubsetbuf(Buffer[i], 4096);
      in[i].open(Dir + '/' + File[i], std::ios_base::binary);
    }
  } else {

    for (uint8_t i = 0; i < 4; ++i)
      in[i].close();
  }
}

void Read_Halos_Buffered (std::array<std::ifstream, 4> &in,
                          std::array<char *, 4> &buffer, uint32_t buff_sz) {
  uint8_t i = 0;
  uint_fast32_t p;

#pragma omp parallel for num_threads(4)
  for (i = 0; i < 4; ++i)
    in[i].read(buffer[i], sizeof(float) * buff_sz);
}
