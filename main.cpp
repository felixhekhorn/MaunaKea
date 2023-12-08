#include "./MaunaKea.hpp"

int main() {
  cuint nl = 5;
  cdbl m2 = pow(172.5, 2);
  cdbl S_h = pow(8e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, 3, 7);
  mk.setHadronicS(S_h);
  mk.setPDF("NNPDF40_nnlo_as_01180", 0);
  // fill the grid
  mk.run();
  printf("res: %e +- %e [pb]\n", mk.intOut.result, mk.intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}
