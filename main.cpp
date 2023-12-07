#include <LHAPDF/LHAPDF.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
#include <vector>

#include "./MaunaKea.hpp"

int main() {
  cuint nl = 5;
  cdbl m2 = pow(172.5, 2);
  cdbl S_h = pow(8e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, 1, 3);
  mk.setHadronicS(S_h);
  mk.setPDF("MSTW2008nnlo68cl", 0);
  // fill the grid
  mk.run();
  printf("res: %e +- %e [pb]\n", mk.intOut.result, mk.intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}
