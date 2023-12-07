#ifndef KERNEL_HPP_
#define KERNEL_HPP_

#include <memory>
#include <vector>

#include "./FO.hpp"
#include "./PineAPPL.hpp"
#include "./config.h"

namespace MaunaKea {
/**
 * @brief Integration kernel
 */
class Kernel : public HepSource::Integrand {
  /** @brief heavy quark mass */
  dbl m2;
  /** @brief number of light flavors */
  uint nl;
  /** @brief active orders by bit */
  uint order_mask;
  /** @brief active flavors by bit */
  uint flavor_mask;
  /** @brief hadronic energy parameter \f$\rho_h = 4m^2/S_h\f$ */
  dbl rho_h;
  /** @brief hadronic energy parameter \f$\beta_h = \sqrt{1-\rho_h}\f$ */
  dbl beta_h;
  /** @brief grid */
  std::unique_ptr<PineAPPL::Grid> grid;
  /** @brief reference PDF */
  std::unique_ptr<LHAPDF::PDF> pdf;
  /** @brief strong coupling constant */
  dbl as;

 public:
  /** @name Order masks */
  ///@{
  /** @brief LO marker */
  static cuint ORDER_LO = 1;
  /** @brief NLO marker */
  static cuint ORDER_NLO = 1 << 1;
  /** @brief NNLO marker */
  static cuint ORDER_NNLO = 1 << 2;
  /** @brief all order marker */
  static cuint ORDER_ALL = ORDER_LO | ORDER_NLO | ORDER_NNLO;
  ///@}

  /** @name Flavor masks */
  ///@{
  /** @brief gluon-gluon marker */
  static cuint FLAVOR_GG = 1;
  /** @brief quark-antiquark marker */
  static cuint FLAVOR_QQBAR = 1 << 1;
  /** @brief all order marker */
  static cuint FLAVOR_ALL = FLAVOR_GG | FLAVOR_QQBAR;
  ///@}

  /**
   * @brief Constructor
   * @param m2 heavy quark mass
   * @param nl number of light flavors
   */
  Kernel(cdbl m2, cuint nl, cuint order_mask, cuint flavor_mask)
      : m2(m2), nl(nl), order_mask(order_mask), flavor_mask(flavor_mask) {}

  /**
   * @brief Set hadronic Mandelstam \f$S_h\f$
   * @param Sh hadronic Mandelstam S
   */
  void setHadronicS(cdbl Sh) {
    this->rho_h = 4. * this->m2 / Sh;
    this->beta_h = sqrt(1. - this->rho_h);
  }

  /**
   * @brief Set reference PDF
   * @param setname PDF set name
   * @param member PDF member
   */
  void setPDF(const str setname, cuint member) {
    this->pdf.reset(LHAPDF::mkPDF(setname, member));
    this->as = this->pdf->alphasQ2(this->m2);
  }

  /**
   * @brief was a PDF set?
   * @return PDF set?
   */
  bool hasPDF() const { return static_cast<bool>(this->pdf); }

  /**
   * @brief Kernel function in Dvegas
   * @param x adapted continuous integration variables
   * @param k discrete integration variables
   * @param vegas_weight integration weight
   * @param aux unadapted continuous integration variables
   * @param f output
   */
  void operator()(const double x[], const int k[], const double& vegas_weight, cdbl aux[], double f[]) {
    (void)k;
    (void)aux;
    cdbl beta_min = 0.;
    cdbl beta_max = this->beta_h;
    cdbl beta = beta_min + (beta_max - beta_min) * x[0];
    cdbl rho = 1 - beta * beta;
    cdbl tau = rho_h / rho;
    cdbl x1_min = tau;
    cdbl x1_max = 1.;
    cdbl x1 = x1_min + (x1_max - x1_min) * x[1];
    cdbl x2 = tau / x1;
    cdbl raw_jac = 2. * beta / rho / x1;
    cdbl jac = raw_jac * (beta_max - beta_min) * (x1_max - x1_min);
    // 1/m2 to get the dimension correct + convert to pb
    cdbl norm = 0.38937966e9 / this->m2;
    cdbl common_weight = norm * jac;
    cdbl mu2 = this->m2;
    // Collect all pieces
    dbl tot = 0.;
    // gluon-gluon channel
    if ((this->flavor_mask & FLAVOR_GG) == FLAVOR_GG) {
      cdbl gg = this->pdf->xfxQ2(21, x1, mu2) * this->pdf->xfxQ2(21, x2, mu2);
      if ((this->order_mask & ORDER_LO) == ORDER_LO) {
        cdbl weight = common_weight * f0gg(rho);
        this->grid->fill(x1, x2, mu2, 0, 0.5, 0, weight * vegas_weight * x1 * x2);
        tot += weight * gg * pow(this->as, 2);
      }
    }
    f[0] = tot;
  }

  /** @brief Initialize grid */
  void initGrid() {
    // create a new luminosity function
    PineAPPL::Lumi lumi;
    lumi.add({PineAPPL::LumiEntry{21, 21, 1.0}});

    // only LO
    std::vector<PineAPPL::Order> orders = {PineAPPL::Order{2, 0, 0, 0}};

    // fully-inclusive cross-section
    std::vector<double> bins = {0.0, 1.0};

    // create the PineAPPL grid
    PineAPPL::KeyVal kv;
    kv.set_double("x_min", rho_h);
    this->grid.reset(new PineAPPL::Grid(lumi, orders, bins, kv));
  }

  /** @brief Optimize grid */
  void optimizeGrid() const { this->grid->optimize(); }

  /** @brief Write grid to disk */
  void writeGrid(const str fp) const { this->grid->write(fp); }

  /** @see HepSource::Integrand::Dvegas_init */
  void Dvegas_init() const { this->grid->scale(0.); }

  /** @see HepSource::Integrand::Dvegas_final */
  void Dvegas_final(cuint iterations) const { this->grid->scale(1. / iterations); }
};
}  // namespace MaunaKea

// if ((p->partonchannel)=="qqbar")
// {
// f = (p->Flux).fluxqqbar(x) * (p->FO).FOqqbar((p->rho)/x) ;
// }
// else if ((p->partonchannel)=="gg")
// {
// f = (p->Flux).fluxgg(x) * (p->FO).FOgg((p->rho)/x) ;
// }
// else if ((p->partonchannel)=="qg")
// {
// f = (p->Flux).fluxgq(x) * (p->FO).FOgq((p->rho)/x) ;
// }
// else if ((p->partonchannel)=="qq")
// {
// f = (p->Flux).fluxqq(x) * (p->FO).FOqq((p->rho)/x) ;
// }
// else if ((p->partonchannel)=="qqprime")
// {
// f = (p->Flux).fluxqqprime(x) * (p->FO).FOqqprime((p->rho)/x) ;
// }
// else if ((p->partonchannel)=="qqbarprime")
// {
// f = (p->Flux).fluxqqbarprime(x) * (p->FO).FOqqbarprime((p->rho)/x) ;
// }
// else//all channels combined
// {
// f = (p->Flux).fluxqqbar(x)    * (p->FO).FOqqbar((p->rho)/x)      +
// (p->Flux).fluxgg(x)         * (p->FO).FOgg((p->rho)/x)         +
// (p->Flux).fluxgq(x)         * (p->FO).FOgq((p->rho)/x)         +
// (p->Flux).fluxqq(x)         * (p->FO).FOqq((p->rho)/x)         +
// (p->Flux).fluxqqprime(x)    * (p->FO).FOqqprime((p->rho)/x)    +
// (p->Flux).fluxqqbarprime(x) * (p->FO).FOqqbarprime((p->rho)/x) ;
// }
// return f/x;
// }

#endif  // KERNEL_HPP_
