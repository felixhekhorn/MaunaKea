#ifndef KERNEL_HPP_
#define KERNEL_HPP_

#include <memory>
#include <vector>

#include "./FO.hpp"
#include "./PineAPPL.hpp"
#include "./config.h"

namespace MaunaKea {

/** @brief Integration kernel */
class Kernel : public HepSource::Integrand {
  /** @brief heavy quark mass */
  dbl m2;
  /** @brief number of light flavors */
  uint nl;
  /** @brief active orders by bit */
  uint order_mask;
  /** @brief active luminosities by bit */
  uint lumi_mask;
  /** @brief hadronic energy parameter \f$\rho_h = 4m^2/S_h\f$ */
  dbl rho_h;
  /** @brief hadronic energy parameter \f$\beta_h = \sqrt{1-\rho_h}\f$ */
  dbl beta_h;
  /** @brief grid */
  std::unique_ptr<PineAPPL::Grid> grid;
  /** @brief reference PDF */
  std::unique_ptr<LHAPDF::PDF> pdf;
  /** @brief factorization scale */
  dbl muF2;
  /** @brief strong coupling constant */
  dbl as;

  /** @brief Integration variables */
  struct IntVars {
    cdbl beta_min = 0.;
    cdbl x1_max = 1.;
    dbl beta_max = dblNaN;
    dbl rho = dblNaN;
    dbl rho_h = dblNaN;
    dbl x1 = dblNaN;
    dbl x2 = dblNaN;
    dbl jac = dblNaN;
    dbl norm = dblNaN;
    dbl common_weight = dblNaN;
    dbl vegas_weight = dblNaN;

    void update(cdbl a0, cdbl a1, cdbl vegas_weight) {
      cdbl beta = this->beta_min + (this->beta_max - this->beta_min) * a0;
      this->rho = 1 - beta * beta;
      cdbl tau = rho_h / rho;
      cdbl x1_min = tau;
      this->x1 = x1_min + (this->x1_max - x1_min) * a1;
      this->x2 = tau / x1;
      cdbl raw_jac = 2. * beta / this->rho / this->x1;
      cdbl jac = raw_jac * (this->beta_max - this->beta_min) * (this->x1_max - x1_min);
      this->common_weight = this->norm * jac;
      this->vegas_weight = vegas_weight;
    }
  };

  /** @brief integration variables */
  IntVars v;

  /** @name Order mapping */
  ///@{
  /** @brief LO position */
  static cuint IDX_ORDER_LO = 0;
  /** @brief NLO position */
  static cuint IDX_ORDER_NLO = 1;
  /** @brief NLO ren. SV position */
  static cuint IDX_ORDER_NLO_R = 2;
  /** @brief NLO fact. SV position */
  static cuint IDX_ORDER_NLO_F = 3;
  ///@}

  /** @name Luminosity mapping */
  ///@{
  /** @brief gluon-gluon position */
  static cuint IDX_LUMI_GG = 0;
  /** @brief quark-antiquark position */
  static cuint IDX_LUMI_QQBAR = 1;
  /** @brief gluon-quark position */
  static cuint IDX_LUMI_GQ = 2;
  ///@}

  /** @brief hide inherited assignment */
  using HepSource::Integrand::operator=;

  /** @brief gluon-gluon flux */
  dbl flux_gg() const {
    cdbl gg = this->pdf->xfxQ2(21, this->v.x1, this->muF2) * this->pdf->xfxQ2(21, this->v.x2, this->muF2);
    return gg;
  }

  /** @brief quark-antiquark flux */
  dbl flux_qqbar() const {
    dbl qqbar = 0.;
    for (uint pid = 1; pid <= this->nl; ++pid) {
      qqbar += 2. * this->pdf->xfxQ2(pid, this->v.x1, this->muF2) * this->pdf->xfxQ2(-pid, this->v.x2, this->muF2);
    }
    return qqbar;
  }

  /** @brief gluon-quark flux */
  dbl flux_gq() const {
    dbl gq = 0.;
    for (uint pid = 1; pid <= this->nl; ++pid) {
      gq += 2. * this->pdf->xfxQ2(21, this->v.x1, this->muF2) * this->pdf->xfxQ2(-pid, this->v.x2, this->muF2);
      gq += 2. * this->pdf->xfxQ2(21, this->v.x1, this->muF2) * this->pdf->xfxQ2(pid, this->v.x2, this->muF2);
    }
    return gq;
  }

  /**
   * @brief Insert a given luminosity into the grid
   * @param flux flux
   * @param idx_lumi grid luminosity index
   * @param m coefficient functions
   */
  dbl fillLumi(cdbl flux, cuint idx_lumi, const FixedOrder::CoeffMap m) const {
    dbl tot = 0.;
    // LO
    if ((this->order_mask & ORDER_LO) == ORDER_LO) {
      cdbl weight = this->v.common_weight * m.f0(this->v.rho, this->nl);
      this->grid->fill(this->v.x1, this->v.x2, this->muF2, IDX_ORDER_LO, 0.5, idx_lumi,
                       weight * this->v.vegas_weight * this->v.x1 * this->v.x2);
      tot += weight * flux * pow(this->as, 2);
    }
    // NLO
    if ((this->order_mask & ORDER_NLO) == ORDER_NLO) {
      if (m.f1) {  // bare
        cdbl weight = this->v.common_weight * m.f1(this->v.rho, this->nl);
        this->grid->fill(this->v.x1, this->v.x2, this->muF2, IDX_ORDER_NLO, 0.5, idx_lumi,
                         weight * this->v.vegas_weight * this->v.x1 * this->v.x2);
        tot += weight * flux * pow(this->as, 3);
      }
      if (m.fbarR1) {  // R SV
        cdbl weight = this->v.common_weight * m.fbarR1(this->v.rho, this->nl);
        this->grid->fill(this->v.x1, this->v.x2, this->muF2, IDX_ORDER_NLO_R, 0.5, idx_lumi,
                         weight * this->v.vegas_weight * this->v.x1 * this->v.x2);
        // tot += weight * flux * pow(this->as, 3);
      }
      if (m.fbarF1) {  // F SV
        cdbl weight = this->v.common_weight * m.fbarF1(this->v.rho, this->nl);
        this->grid->fill(this->v.x1, this->v.x2, this->muF2, IDX_ORDER_NLO_F, 0.5, idx_lumi,
                         weight * this->v.vegas_weight * this->v.x1 * this->v.x2);
        // tot += weight * flux * pow(this->as, 3);
      }
    }
    return tot;
  }

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

  /** @name Luminosity masks */
  ///@{
  /** @brief gluon-gluon marker */
  static cuint LUMI_GG = 1;
  /** @brief quark-antiquark marker */
  static cuint LUMI_QQBAR = 1 << 1;
  /** @brief gluon-quark marker */
  static cuint LUMI_GQ = 1 << 2;
  /** @brief all order marker */
  static cuint LUMI_ALL = LUMI_GG | LUMI_QQBAR | LUMI_GQ;
  ///@}

  /**
   * @brief Constructor
   * @param m2 heavy quark mass
   * @param nl number of light flavors
   * @param order_mask active orders
   * @param lumi_mask active luminosities
   */
  Kernel(cdbl m2, cuint nl, cuint order_mask, cuint lumi_mask)
      : m2(m2), nl(nl), order_mask(order_mask), lumi_mask(lumi_mask) {
    // 1/m2 to get the dimension correct and convert to pb
    this->v.norm = 0.38937966e9 / this->m2;
    this->muF2 = this->m2;
  }

  /**
   * @brief Set hadronic Mandelstam \f$S_h\f$
   * @param Sh hadronic Mandelstam S
   */
  void setHadronicS(cdbl Sh) {
    this->rho_h = 4. * this->m2 / Sh;
    this->beta_h = sqrt(1. - this->rho_h);
    this->v.beta_max = this->beta_h;
    this->v.rho_h = this->rho_h;
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
  void operator()(const double x[], const int k[], const double& vegas_weight, const double aux[], double f[]) {
    (void)k;
    (void)aux;
    this->v.update(x[0], x[1], vegas_weight);
    // Collect all pieces
    dbl tot = 0.;
    // gluon-gluon channel
    if ((this->lumi_mask & LUMI_GG) == LUMI_GG) tot += this->fillLumi(this->flux_gg(), IDX_LUMI_GG, FixedOrder::gg);
    // quark-antiquark channel
    if ((this->lumi_mask & LUMI_QQBAR) == LUMI_QQBAR)
      tot += this->fillLumi(this->flux_qqbar(), IDX_LUMI_QQBAR, FixedOrder::qqbar);
    // gluon-quark channel
    if ((this->lumi_mask & LUMI_GQ) == LUMI_GQ) tot += this->fillLumi(this->flux_gq(), IDX_LUMI_GQ, FixedOrder::gq);
    f[0] = tot;
  }

  /** @brief Initialize grid */
  void initGrid() {
    // create a new luminosity function
    PineAPPL::Lumi lumi;
    // gluon-gluon channel
    lumi.add({PineAPPL::LumiEntry{21, 21, 1.0}});
    // quark-antiquark channel
    {
      std::vector<PineAPPL::LumiEntry> qqbar;
      for (uint pid = 1; pid <= this->nl; ++pid) qqbar.push_back({static_cast<int>(pid), static_cast<int>(-pid), 2.0});
      lumi.add(qqbar);
    }
    // gluon-quark channel
    {
      std::vector<PineAPPL::LumiEntry> gq;
      for (uint pid = 1; pid <= this->nl; ++pid) {
        gq.push_back({21, static_cast<int>(-pid), 2.0});
        gq.push_back({21, static_cast<int>(pid), 2.0});
      }
      lumi.add(gq);
    }

    std::vector<PineAPPL::Order> orders = {// LO
                                           PineAPPL::Order{2, 0, 0, 0},
                                           // NLO
                                           PineAPPL::Order{3, 0, 0, 0}, PineAPPL::Order{3, 0, 1, 0},
                                           PineAPPL::Order{3, 0, 0, 1}};

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
