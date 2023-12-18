#include "./MaunaKea.hpp"

int toppp() {
  cuint nl = 5;
  cdbl m2 = pow(172.5, 2);
  cdbl S_h = pow(8e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_ALL, MaunaKea::Kernel::LUMI_ALL);
  mk.intCfg.calls = 50000;
  mk.setHadronicS(S_h);
  mk.setPDF("NNPDF40_nnlo_as_01180", 0);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}

int dEnterria_charm() {
  cuint nl = 3;
  cdbl m2 = pow(1.67, 2);
  cdbl S_h = pow(14e2, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_LO, MaunaKea::Kernel::LUMI_GG);
  mk.intCfg.calls = 50000;
  mk.setHadronicS(S_h);
  mk.setPDF("NNPDF40_nnlo_as_01180", 0);
  // mk.setScaleRatios(2.,2.);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}

int dEnterria_bottom() {
  cuint nl = 4;
  cdbl m2 = pow(4.66, 2);
  cdbl S_h = pow(14e2, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_LO, MaunaKea::Kernel::LUMI_GG);
  mk.intCfg.calls = 50000;
  mk.setHadronicS(S_h);
  mk.setPDF("NNPDF40_nnlo_as_01180", 0);
  // mk.setScaleRatios(2.,2.);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea.pineappl.lz4");
  return EXIT_SUCCESS;
}

int main() {
  // return toppp();
  // return dEnterria_charm();
  return dEnterria_bottom();
}
