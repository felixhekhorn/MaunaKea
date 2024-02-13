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

int toppp() {
  cuint nl = 5;
  cdbl m2 = pow(172.5, 2);
  cdbl S_h = pow(7e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::ORDER_LO, MaunaKea::LUMI_GG);
  // mk.intCfg.calls = 5000;
  mk.intCfg.verbosity = 3;
  mk.setHadronicS(S_h);
  mk.setPDF("NNPDF40_nnlo_as_01180", 0);
  // mk.setCentralScaleRatio(2.);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.getIntegrationOutput();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}

int main() {
  // return Hathor();
  return toppp();
}