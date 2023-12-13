#include "./MaunaKea.hpp"

int main() {
  cuint nl = 5;
  cdbl m2 = pow(172.5, 2);
  cdbl S_h = pow(8e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_NNLO, MaunaKea::Kernel::LUMI_QQPRIME);
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
