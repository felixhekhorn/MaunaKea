
#ifndef Integration_HPP_
#define Integration_HPP_

#include <dvegas/dvegas.h>

#include "config.h"

/**
 * @class IntegrationConfig
 * @brief configurates a single integration
 */
struct IntegrationConfig {
    
/** @name common variables */
///@{
    
/** @brief level of output */
    int verbosity = 0;

/** @brief calls */
    size_t calls = 0;
///@}
    
/** @name Monte Carlo variables */
///@{
/** @brief calls for warmup */
    size_t MC_warmupCalls = 0;
    
/** @brief iterations */
    uint MC_iterations = 5;
    
/** @brief iterations during warmup */
    uint MC_warmupIterations = 5;
    
/** @brief iterate until |chi2-1| < 0.5? */
    bool MC_adaptChi2 = true;
///@}

/** @name variables for Dvegas */
///@{
/** @brief number of bins */
    uint Dvegas_bins = 250;
///@}

};

/**
 * @class IntegrationOutput
 * @brief stores meta data for a single integration
 */
struct IntegrationOutput {
    
/** @name common variables */
///@{
/** @brief result */
    dbl result = 0;
/** @brief (absolute) error */
    dbl error = 0;
///@}
    
/** @name Monte Carlo variables */
///@{
/** @brief chi^2 */
    dbl MC_chi2 = 0;
/** @brief number of iteration to converge chi^2 */
    uint MC_chi2inter = 0;
///@}
    
/**
 * @brief constructor
 * @param result result
 * @param error (absolute) error
 * @param MC_chi2 chi^2
 * @param MC_chi2inter number of iteration to converge chi^2
 */
    IntegrationOutput(cdbl result = 0, cdbl error = 0, cdbl MC_chi2 = 0, cuint MC_chi2inter = 0) : 
        result(result), error(error), MC_chi2(MC_chi2), MC_chi2inter(MC_chi2inter) {};

};

/**
 * @brief integrates the kernel in 2 dimension
 * @param K kernel
 * @param cfg config
 * @param out output of meta data
 * @return \f$\int\limits_0^1 f(a_1,a_2)\,da_1da_2\f$
 */
template <class IntKerT> dbl integrate2D(IntKerT* K, const IntegrationConfig& cfg, IntegrationOutput* out) {
    cuint dim = 2;
    HepSource::Dvegas dv(dim,cfg.Dvegas_bins,1,*K);
    /** @todo activate correlation between z and x? -> Dvegas dv(dim,cfg.Dvegas_bins,2,{},0,1,F); */
    double res, err;
    /* clear histograms */
    K->Dvegas_init();
    /* warm-up */
    /* catch zero kernel */
    try{ HepSource::VEGAS(dv,cfg.MC_warmupCalls,cfg.MC_warmupIterations,0,cfg.verbosity - 3); }
    catch(domain_error& e){ out->result = 0; out->error = 0; out->MC_chi2 = 0; out->MC_chi2inter = 0; return 0.; }
    HepSource::IntegrandEstimate e = dv.stats(0);
    res = e.integral();
    uint guard = 0;
    out->result = dblNaN;
    out->error = dblNaN;
    out->MC_chi2 = dblNaN;
    out->MC_chi2inter = 0;
    /* run */
    if (cfg.MC_adaptChi2) { /* adapt chi */
        do {
            if (!isfinite(res)) return res;
            K->Dvegas_init();
            HepSource::VEGAS(dv,cfg.calls,cfg.MC_iterations,1,cfg.verbosity - 2);
            e = dv.stats(0);
            res = e.integral();
            err = e.standardDeviation();
            if (cfg.verbosity > 1)
                printf("[INFO] int%dD(Dvegas): [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n",dim,guard,res,err,fabs(err/res*1e2),e.chiSquarePerIteration());
        } while (fabs (e.chiSquarePerIteration() - 1.0) > 0.5 && ++guard < 15);
    } else { /* simple run */
        K->Dvegas_init();
        HepSource::VEGAS(dv,cfg.calls,cfg.MC_iterations,1,cfg.verbosity - 2);
        e = dv.stats(0);
        res = e.integral();
        err = e.standardDeviation();
    }
    K->Dvegas_final(cfg.MC_iterations);
    if (1 == cfg.verbosity)
        printf("[INFO] int%dD(Dvegas): [%d] % e ± %.3e (%.3f%%) chi2/it: %.3f\n",dim,guard,res,err,fabs(err/res*1e2),e.chiSquarePerIteration());
    out->result = res;
    out->error = err;
    out->MC_chi2 = e.chiSquarePerIteration();
    out->MC_chi2inter = guard;
    return res;
}

#endif // Integration_HPP_
