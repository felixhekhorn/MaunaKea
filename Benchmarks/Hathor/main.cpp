#include "./MaunaKea.hpp"

int Hathor() {
  cuint nl = 3;
  cdbl m2 = pow(1.51, 2);
  cdbl S_h = pow(8e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::ORDER_LO | MaunaKea::ORDER_NLO, MaunaKea::LUMI_ALL);
  // mk.intCfg.calls = 100000;
  mk.intCfg.verbosity = 3;
  mk.setHadronicS(S_h);
  mk.setPDF("NNPDF40_nlo_pch_as_01180_nf_3", 0);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.getIntegrationOutput();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}

int main() { return Hathor(); }
