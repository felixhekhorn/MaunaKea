#include "./MaunaKea.hpp"

int toppp() {
  cuint nl = 5;
  cdbl m2 = pow(172.5, 2);
  cdbl S_h = pow(8e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_LO, MaunaKea::Kernel::LUMI_GG);
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

int dEnterria_table1_charm() {
  cuint nl = 3;
  cdbl m2 = pow(1.67, 2);
  // cdbl S_h = pow(14e3, 2);
  cdbl S_h = pow(100e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_ALL, MaunaKea::Kernel::LUMI_ALL);
  mk.intCfg.calls = 50000;
  mk.setHadronicS(S_h);
  mk.setPDF("ABMP16_3_nnlo", 0);
  mk.setScaleRatios(2., 2.);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea-ccbar.pineappl.lz4");
  return EXIT_SUCCESS;
}

int dEnterria_fig1_charm() {
  cuint nl = 3;
  cdbl m2 = pow(1.67, 2);
  cdbl Sh_min = pow(35., 2);
  cdbl Sh_max = pow(100e3, 2);
  cuint ndata = 10;
  for (uint j = 0; j < ndata; ++j) {
    cdbl logS_h = log(Sh_min) + (log(Sh_max) - log(Sh_min)) * j / (ndata - 1);
    cdbl S_h = exp(logS_h);
    printf("j = %d, sqrt(S) = %e", j, S_h);
    // init object
    MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_ALL, MaunaKea::Kernel::LUMI_ALL);
    mk.intCfg.calls = 50000;
    mk.setHadronicS(S_h);
    mk.setPDF("ABMP16_3_nnlo", 0);
    mk.setScaleRatios(2., 2.);
    // fill the grid
    mk.run();
    const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
    printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
    // save
#define kPathSize 100
    char buffer[kPathSize];
    snprintf(buffer, kPathSize, "MaunaKea-fig1-ccbar-%d.pineappl.lz4", j);
    mk.write(buffer);
  }
  return EXIT_SUCCESS;
}

int dEnterria_table1_bottom() {
  cuint nl = 4;
  cdbl m2 = pow(4.66, 2);
  // cdbl S_h = pow(14e3, 2);
  cdbl S_h = pow(100e3, 2);
  // init object
  MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_ALL, MaunaKea::Kernel::LUMI_ALL);
  mk.intCfg.calls = 50000;
  mk.setHadronicS(S_h);
  mk.setPDF("ABMP16_4_nnlo", 0);
  mk.setScaleRatios(2., 2.);
  // fill the grid
  mk.run();
  const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
  printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
  // save
  mk.write("MaunaKea-bbbar.pineappl.lz4");
  return EXIT_SUCCESS;
}

int dEnterria_fig1_bottom() {
  cuint nl = 4;
  cdbl m2 = pow(4.66, 2);
  cdbl Sh_min = pow(35., 2);
  cdbl Sh_max = pow(100e3, 2);
  cuint ndata = 10;
  for (uint j = 0; j < ndata; ++j) {
    cdbl logS_h = log(Sh_min) + (log(Sh_max) - log(Sh_min)) * j / (ndata - 1);
    cdbl S_h = exp(logS_h);
    printf("j = %d, sqrt(S) = %e", j, S_h);
    // init object
    MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_ALL, MaunaKea::Kernel::LUMI_ALL);
    mk.intCfg.calls = 50000;
    mk.setHadronicS(S_h);
    mk.setPDF("ABMP16_4_nnlo", 0);
    mk.setScaleRatios(2., 2.);
    // fill the grid
    mk.run();
    const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
    printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
    // save
#define kPathSize 100
    char buffer[kPathSize];
    snprintf(buffer, kPathSize, "MaunaKea-fig1-bbbar-%d.pineappl.lz4", j);
    mk.write(buffer);
  }
  return EXIT_SUCCESS;
}

int main() {
  // return toppp();
  // return dEnterria_table1_charm();
  // return dEnterria_fig1_charm();
  // return dEnterria_table1_bottom();
  return dEnterria_fig1_bottom();
}
