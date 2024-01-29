#include "../MaunaKea.hpp"

int fig1() {
  cuint nl = 3;
  cdbl m2 = pow(1.51, 2);
  cdbl Sh_min = pow(20., 2);
  cdbl Sh_max = pow(1e4, 2);
  cuint ndata = 25;
  for (uint j = 0; j < ndata; ++j) {
    cdbl logS_h = log(Sh_min) + (log(Sh_max) - log(Sh_min)) * j / (ndata - 1);
    cdbl S_h = exp(logS_h);
    printf("j = %d, sqrt(S) = %e\n", j, pow(S_h, 0.5));
    // init object
    MaunaKea::MaunaKea mk(m2, nl, MaunaKea::Kernel::ORDER_ALL, MaunaKea::Kernel::LUMI_ALL);
    mk.intCfg.calls = 50000;
    mk.setHadronicS(S_h);
    mk.setPDF("NNPDF40_nlo_pch_as_01180_nf_3", 0);
    mk.setGridCentralScaleRatio(2.);
    // fill the grid
    mk.run();
    const MaunaKea::IntegrationOutput intOut = mk.get_integration_output();
    printf("sigma_tot = %e +- %e [pb]\n", intOut.result, intOut.error);
    // save
#define kPathSize 100
    char buffer[kPathSize];
    snprintf(buffer, kPathSize, "MaunaKea-ccbar-fig1-%d.pineappl.lz4", j);
    mk.write(buffer);
  }
  return EXIT_SUCCESS;
}

int main() { return fig1(); }
