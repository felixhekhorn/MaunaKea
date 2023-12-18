
#ifndef INTEGRATION_HPP_
#define INTEGRATION_HPP_

#include <dvegas/dvegas.h>

#include "./config.h"

namespace MaunaKea {

/** @brief Integration parameter configuration */
struct IntegrationConfig {
  /** @brief level of output */
  int verbosity = 0;
  /** @brief calls */
  size_t calls = 0;
  /** @brief calls for warmup */
  size_t warmupCalls = 0;
  /** @brief iterations */
  uint iterations = 5;
  /** @brief iterations during warmup */
  uint warmupIterations = 5;
  /** @brief iterate until |chi2-1| < 0.5? */
  bool adaptChi2 = true;
  /** @brief number of bins */
  uint bins = 250;

  /**
   * @brief Dump current state to string
   * @return string representation
   */
  str toString() const {
#define kIntegrationConfigStrSize = 200
    char buffer[kIntegrationConfigStrSize];
    snprintf(buffer, kIntegrationConfigStrSize,
             "verbosity: %d\ncalls: %lu\nwarmupCalls: %lu\niterations: %u\nwarmupIterations: %u\nadaptChi2:%u\nbins:%u",
             this->verbosity, this->calls, this->warmupCalls, this->iterations, this->warmupIterations, this->adaptChi2,
             this->bins);
    return str(buffer);
  }
};

/** @brief Integration result */
struct IntegrationOutput {
  /** @brief result */
  dbl result = 0;
  /** @brief (absolute) error */
  dbl error = 0;
  /** @brief chi^2 */
  dbl chi2 = 0;
  /** @brief number of iteration to converge chi^2 */
  uint chi2iter = 0;

  /**
   * @brief constructor
   * @param result result
   * @param error (absolute) error
   * @param chi2 chi^2
   * @param chi2iter number of iteration to converge chi^2
   */
  explicit IntegrationOutput(cdbl result = 0, cdbl error = 0, cdbl chi2 = 0, cuint chi2iter = 0)
      : result(result), error(error), chi2(chi2), chi2iter(chi2iter) {}

  /**
   * @brief Dump current state to string
   * @return string representation
   */
  str toString() const {
#define kIntegrationOutputStrSize 200
    char buffer[kIntegrationOutputStrSize];
    snprintf(buffer, kIntegrationOutputStrSize, "result: %e\nerror: %e\nchi2: %e\nchi2iter: %u", this->result,
             this->error, this->chi2, this->chi2iter);
    return str(buffer);
  }
};

/**
 * @brief integrates the kernel in 2 dimension
 * @param K kernel
 * @param cfg configuration
 * @param out output
 * @return \f$\int\limits_0^1 K(a_1,a_2)\,da_1da_2\f$
 */
template <class IntKerT>
dbl integrate2D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out) {
  cuint dim = 2;
  HepSource::Dvegas dv(dim, cfg.bins, 1, *K);
  double res, err;
  // clear histograms
  K->Dvegas_init();
  // warm-up
  // catch zero kernel
  try {
    HepSource::VEGAS(dv, cfg.warmupCalls, cfg.warmupIterations, 0, cfg.verbosity - 3);
  } catch (std::domain_error& e) {
    out->result = 0;
    out->error = 0;
    out->chi2 = 0;
    out->chi2iter = 0;
    return 0.;
  }
  HepSource::IntegrandEstimate e = dv.stats(0);
  res = e.integral();
  uint guard = 0;
  out->result = dblNaN;
  out->error = dblNaN;
  out->chi2 = dblNaN;
  out->chi2iter = 0;
  // run
  if (cfg.adaptChi2) {  // adapt chi
    do {
      if (!std::isfinite(res)) return res;
      K->Dvegas_init();
      HepSource::VEGAS(dv, cfg.calls, cfg.iterations, 1, cfg.verbosity - 2);
      e = dv.stats(0);
      res = e.integral();
      err = e.standardDeviation();
      if (cfg.verbosity > 1)
        printf("[INFO] int%dD(Dvegas): [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n", dim, guard, res, err,
               fabs(err / res * 1e2), e.chiSquarePerIteration());
    } while (fabs(e.chiSquarePerIteration() - 1.0) > 0.5 && ++guard < 15);
  } else {  // simple run
    K->Dvegas_init();
    HepSource::VEGAS(dv, cfg.calls, cfg.iterations, 1, cfg.verbosity - 2);
    e = dv.stats(0);
    res = e.integral();
    err = e.standardDeviation();
  }
  // finish
  K->Dvegas_final(cfg.iterations);
  if (1 == cfg.verbosity)
    printf("[INFO] int%dD(Dvegas): [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n", dim, guard, res, err,
           fabs(err / res * 1e2), e.chiSquarePerIteration());
  out->result = res;
  out->error = err;
  out->chi2 = e.chiSquarePerIteration();
  out->chi2iter = guard;
  return res;
}
}  // namespace MaunaKea

#endif  // INTEGRATION_HPP_
