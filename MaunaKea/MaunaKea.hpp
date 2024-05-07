#ifndef MAUNAKEA_MAUNAKEA_HPP_
#define MAUNAKEA_MAUNAKEA_HPP_

#include "./Integration.hpp"
#include "./Kernel.hpp"

/** @brief Hadroproduction of heavy quark flavors */
namespace MaunaKea {

/** @brief Main application class */
class MaunaKea {
  /** @brief Integration kernel */
  Kernel k;
  /** @brief Integration output */
  IntegrationOutput intOut;

 public:
  /** @brief Integration configuration */
  IntegrationConfig intCfg;

  /**
   * @brief Constructor
   * @param m2 heavy quark mass
   * @param nl number of light flavors
   * @param order_mask active orders
   * @param lumi_mask active luminosities
   */
  MaunaKea(cdbl m2, cuint nl, cuint order_mask, cuint lumi_mask) : k(m2, nl, order_mask, lumi_mask) {
    this->intCfg.warmupCalls = 1000;
    this->intCfg.calls = 20000;
  }

  /**
   * @brief Set hadronic Mandelstam \f$S_h\f$
   * @param Sh hadronic Mandelstam S
   */
  void setHadronicS(cdbl Sh) { this->k.setHadronicS(Sh); }

  /**
   * @brief Set renormalization and factorization scale ratios \f$\xi_{R/F} = \mu_{R/F}/m\f$
   * @param xiR (linear) renormalization scale ratio \f$\xi_R = \mu_R/m\f$
   * @param xiF (linear) factorization scale ratio \f$\xi_F = \mu_F/m\f$
   */
  void setScaleRatios(cdbl xiR, cdbl xiF) { this->k.setScaleRatios(xiR, xiF); }

  /**
   * @brief Set grid central scale ratio \f$\xi = \mu/m\f$
   * @param xi (linear) central scale ratio \f$\xi = \mu/m\f$
   */
  void setCentralScaleRatio(cdbl xi) { this->k.setCentralScaleRatio(xi); }

  /**
   * @brief Set reference PDF
   * @param setname PDF set name
   * @param member PDF member
   */
  void setPDF(const str& setname, cuint member) { this->k.setPDF(setname, member); }

  /**
   * @brief Set reference PDF
   * @param setname PDF set name and member separated by /
   */
  void setPDF(const str& setname_nmem) { this->k.setPDF(setname_nmem); }

  /** @brief Run calculation */
  void run() {
    // Checks
    if (!this->k.isPhysical()) throw std::domain_error("No physical domain (too little energy)!");
    if (!this->k.hasPDF()) throw std::domain_error("No PDF was set so far!");
    // Init
    this->k.initGrid();
    // do it!
    IntegrationOutput out;
    integrate2D(&this->k, this->intCfg, &out);
    this->intOut = out;
    // post-process
    this->k.optimizeGrid();
    this->k.addScopedMetadata("IntegrationConfig", this->intCfg.toString());
    this->k.addLocalMetadata();
    this->k.addScopedMetadata("IntegrationOutput", out.toString());
  }

  /**
   * @brief Write grid to disk
   * @param fp file path
   */
  void write(const str& fp) const { this->k.writeGrid(fp); }

  /**
   * @brief Copy of current integration output
   * @return Integration output
   */
  IntegrationOutput getIntegrationOutput() const { return IntegrationOutput(this->intOut); }
};
}  // namespace MaunaKea
#endif  // MAUNAKEA_MAUNAKEA_HPP_
