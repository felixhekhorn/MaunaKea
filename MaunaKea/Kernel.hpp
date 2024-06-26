#ifndef MAUNAKEA_KERNEL_HPP_
#define MAUNAKEA_KERNEL_HPP_

#include <memory>
#include <vector>

#include <PineAPPL.hpp>

#include "./FO.hpp"
#include "./config.h"

namespace MaunaKea {

/** @name Order masks */
///@{
/** @brief LO marker */
const cuint ORDER_LO = 1;
/** @brief NLO marker */
const cuint ORDER_NLO = 1 << 1;
/** @brief NNLO marker */
const cuint ORDER_NNLO = 1 << 2;
/** @brief all order marker */
const cuint ORDER_ALL = ORDER_LO | ORDER_NLO | ORDER_NNLO;
///@}

/** @name Luminosity masks */
///@{
/** @brief gluon-gluon marker */
const cuint LUMI_GG = 1;
/** @brief quark-antiquark marker */
const cuint LUMI_QQBAR = 1 << 1;
/** @brief gluon-quark marker */
const cuint LUMI_GQ = 1 << 2;
/** @brief quark-quark marker */
const cuint LUMI_QQ = 1 << 3;
/** @brief quark-antiquark' marker */
const cuint LUMI_QQBARPRIME = 1 << 4;
/** @brief quark-quark' marker */
const cuint LUMI_QQPRIME = 1 << 5;
/** @brief all order marker */
const cuint LUMI_ALL = LUMI_GG | LUMI_QQBAR | LUMI_GQ | LUMI_QQ | LUMI_QQBARPRIME | LUMI_QQPRIME;
///@}

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
  /** @brief hadronic c.o.m. energy  \f$S_h\f$ */
  dbl S_h;
  /** @brief hadronic energy parameter \f$\rho_h = 4m^2/S_h\f$ */
  dbl rho_h;
  /** @brief hadronic energy parameter \f$\beta_h = \sqrt{1-\rho_h}\f$ */
  dbl beta_h;
  /** @brief grid */
  std::unique_ptr<PineAPPL::Grid> grid;
  /** @brief reference PDF */
  std::unique_ptr<LHAPDF::PDF> pdf;
  /** @brief grid central scale ratio \f$\xi = \mu/m\f$ */
  dbl xi = 1.;
  /** @brief grid central scale log \f$L_C=2*\log(\xi)\f$ */
  dbl logC = 0.;
  /** @brief grid central scale */
  dbl mu2;
  /** @brief renormalization scale ratio \f$\xi_R = \mu_R/m\f$ */
  dbl xiR = 1.;
  /** @brief renormalization scale log \f$L_R = 2*\log(\xi_R)\f$ */
  dbl logR = 0.;
  /** @brief strong coupling constant */
  dbl as;
  /** @brief factorization scale ratio \f$\xi_F = \mu_F/m\f$ */
  dbl xiF = 1.;
  /** @brief factorization scale log \f$L_F = 2*\log(\xi_F)\f$ */
  dbl logF = 0.;
  /** @brief factorization scale */
  dbl muF2;

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
  /** @brief NLO R SV position */
  static cuint IDX_ORDER_NLO_R = 2;
  /** @brief NLO F SV position */
  static cuint IDX_ORDER_NLO_F = 3;
  /** @brief NNLO position */
  static cuint IDX_ORDER_NNLO = 4;
  /** @brief NNLO R SV position */
  static cuint IDX_ORDER_NNLO_R = 5;
  /** @brief NNLO F SV position */
  static cuint IDX_ORDER_NNLO_F = 6;
  /** @brief NNLO RR SV position */
  static cuint IDX_ORDER_NNLO_RR = 7;
  /** @brief NNLO RF SV position */
  static cuint IDX_ORDER_NNLO_RF = 8;
  /** @brief NNLO FF SV position */
  static cuint IDX_ORDER_NNLO_FF = 9;
  ///@}

  /** @name Luminosity mapping */
  ///@{
  /** @brief gluon-gluon position */
  static cuint IDX_LUMI_GG = 0;
  /** @brief quark-antiquark position */
  static cuint IDX_LUMI_QQBAR = 1;
  /** @brief gluon-quark position */
  static cuint IDX_LUMI_GQ = 2;
  /** @brief quark-quark position */
  static cuint IDX_LUMI_QQ = 3;
  /** @brief quark-antiquark' position */
  static cuint IDX_LUMI_QQBARPRIME = 4;
  /** @brief quark-quark' position */
  static cuint IDX_LUMI_QQPRIME = 5;
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

  /** @brief quark-quark flux */
  dbl flux_qq() const {
    dbl qq = 0.;
    for (uint pid = 1; pid <= this->nl; ++pid) {
      qq += this->pdf->xfxQ2(pid, this->v.x1, this->muF2) * this->pdf->xfxQ2(pid, this->v.x2, this->muF2);
      qq += this->pdf->xfxQ2(-pid, this->v.x1, this->muF2) * this->pdf->xfxQ2(-pid, this->v.x2, this->muF2);
    }
    return qq;
  }

  /** @brief quark-antiquark' flux */
  dbl flux_qqbarprime() const {
    dbl qqbarprime = 0.;
    for (uint pid1 = 1; pid1 <= this->nl; ++pid1) {
      for (uint pid2 = 1; pid2 <= this->nl; ++pid2) {
        if (pid1 == pid2) continue;
        qqbarprime += this->pdf->xfxQ2(pid1, this->v.x1, this->muF2) * this->pdf->xfxQ2(-pid2, this->v.x2, this->muF2);
        qqbarprime += this->pdf->xfxQ2(-pid1, this->v.x1, this->muF2) * this->pdf->xfxQ2(pid2, this->v.x2, this->muF2);
      }
    }
    return qqbarprime;
  }

  /** @brief quark-quark' flux */
  dbl flux_qqprime() const {
    dbl qqprime = 0.;
    for (uint pid1 = 1; pid1 <= this->nl; ++pid1) {
      for (uint pid2 = 1; pid2 <= this->nl; ++pid2) {
        if (pid1 == pid2) continue;
        qqprime += this->pdf->xfxQ2(pid1, this->v.x1, this->muF2) * this->pdf->xfxQ2(pid2, this->v.x2, this->muF2);
        qqprime += this->pdf->xfxQ2(-pid1, this->v.x1, this->muF2) * this->pdf->xfxQ2(-pid2, this->v.x2, this->muF2);
      }
    }
    return qqprime;
  }

  /**
   * @brief Compute grid contribution and integral contribution
   * @param flux flux
   * @param k PTO
   * @param idx_order grid order index
   * @param powR renormalization log power
   * @param powF factorization log power
   * @param idx_lumi grid luminosity index
   * @param weight coefficient function
   * @return contribution to integral
   */
  dbl fill(cdbl flux, cuint k, cuint idx_order, cuint powR, cuint powF, cuint idx_lumi, cdbl weight) const {
    cdbl normed_weight = this->v.common_weight * weight;
    cdbl grid_weight = normed_weight * this->v.vegas_weight * this->v.x1 * this->v.x2;
    this->grid->fill(this->v.x1, this->v.x2, this->mu2, idx_order, 0.5, idx_lumi, grid_weight);
    return normed_weight * flux * pow(this->as, 2 + k) * pow(this->logR, powR) * pow(this->logF, powF);
  }

  /**
   * @brief Insert a given luminosity into the grid
   * @param flux flux
   * @param idx_lumi grid luminosity index
   * @param m coefficient functions
   * @return contribution to integral
   */
  dbl fillLumi(cdbl flux, cuint idx_lumi, const FixedOrder::CoeffMap& m) const {
    dbl tot = 0.;
    // LO
    if ((this->order_mask & ORDER_LO) == ORDER_LO) {
      if (m.f0) {
        cdbl w0 = m.f0(this->v.rho, this->nl);
        tot += this->fill(flux, 0, this->IDX_ORDER_LO, 0, 0, idx_lumi, w0);
      }
    }
    // NLO
    if ((this->order_mask & ORDER_NLO) == ORDER_NLO) {
      if (m.f1) {
        cdbl w1 = m.f1(this->v.rho, this->nl);
        tot += this->fill(flux, 1, this->IDX_ORDER_NLO, 0, 0, idx_lumi, w1);
      }
      // SV
      if (m.fbarR1) {
        cdbl wR1 = m.fbarR1(this->v.rho, this->nl);
        tot += this->fill(flux, 1, this->IDX_ORDER_NLO_R, 1, 0, idx_lumi, wR1);
        tot += this->fill(flux, 1, this->IDX_ORDER_NLO, 0, 0, idx_lumi, wR1 * this->logC);
      }
      if (m.fbarF1) {
        cdbl wF1 = m.fbarF1(this->v.rho, this->nl);
        tot += this->fill(flux, 1, this->IDX_ORDER_NLO_F, 0, 1, idx_lumi, wF1);
        tot += this->fill(flux, 1, this->IDX_ORDER_NLO, 0, 0, idx_lumi, wF1 * this->logC);
      }
    }
    // NNLO
    if ((this->order_mask & ORDER_NNLO) == ORDER_NNLO) {
      if (m.f2) {
        cdbl w2 = m.f2(this->v.rho, this->nl);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO, 0, 0, idx_lumi, w2);
      }
      // SV^1
      if (m.fbarR2) {
        cdbl wR2 = m.fbarR2(this->v.rho, this->nl);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_R, 1, 0, idx_lumi, wR2);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO, 0, 0, idx_lumi, wR2 * this->logC);
      }
      if (m.fbarF2) {
        cdbl wF2 = m.fbarF2(this->v.rho, this->nl);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_F, 0, 1, idx_lumi, wF2);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO, 0, 0, idx_lumi, wF2 * this->logC);
      }
      // SV^2
      if (m.fbarRR2) {
        cdbl wRR2 = m.fbarRR2(this->v.rho, this->nl);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_RR, 2, 0, idx_lumi, wRR2);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO, 0, 0, idx_lumi, wRR2 * pow(this->logC, 2));
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_R, 1, 0, idx_lumi, wRR2 * 2. * this->logC);
      }
      if (m.fbarRF2) {
        cdbl wRF2 = m.fbarRF2(this->v.rho, this->nl);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_RF, 2, 0, idx_lumi, wRF2);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO, 0, 0, idx_lumi, wRF2 * pow(this->logC, 2));
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_R, 1, 0, idx_lumi, wRF2 * this->logC);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_F, 0, 1, idx_lumi, wRF2 * this->logC);
      }
      if (m.fbarFF2) {
        cdbl wFF2 = m.fbarFF2(this->v.rho, this->nl);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_FF, 0, 2, idx_lumi, wFF2);
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO, 0, 0, idx_lumi, wFF2 * pow(this->logC, 2));
        tot += this->fill(flux, 2, this->IDX_ORDER_NNLO_F, 0, 1, idx_lumi, wFF2 * 2. * this->logC);
      }
    }
    return tot;
  }

  /** @brief Cache alpha_s at the renormalization scale */
  void setAlphaS() {
    cdbl muR2 = this->xiR * this->xiR * this->mu2;
    this->as = this->pdf->alphasQ2(muR2);
  }

  /**
   * @brief Set reference PDF (and alpha_s)
   * @param pdf PDF
   */
  void resetPDF(LHAPDF::PDF* pdf) {
    this->pdf.reset(pdf);
    this->setAlphaS();
  }

 public:
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
    this->setCentralScaleRatio();
    this->setScaleRatios();
  }

  /**
   * @brief Set hadronic Mandelstam \f$S_h\f$
   * @param Sh hadronic Mandelstam S
   */
  void setHadronicS(cdbl Sh) {
    this->S_h = Sh;
    this->rho_h = 4. * this->m2 / Sh;
    this->beta_h = sqrt(1. - this->rho_h);
    this->v.beta_max = this->beta_h;
    this->v.rho_h = this->rho_h;
  }

  /**
   * @brief Set renormalization and factorization scale ratios \f$\xi_{R/F} = \mu_{R/F}/m\f$
   * @param xiR (linear) renormalization scale ratio \f$\xi_R = \mu_R/m\f$
   * @param xiF (linear) factorization scale ratio \f$\xi_F = \mu_F/m\f$
   */
  void setScaleRatios(cdbl xiR = 1., cdbl xiF = 1.) {
    this->xiR = xiR;
    this->logR = 2. * log(xiR);
    // recalculate as
    if (this->hasPDF()) this->setAlphaS();
    this->xiF = xiF;
    this->muF2 = xiF * xiF * this->mu2;
    this->logF = 2. * log(xiF);
  }

  /**
   * @brief Set grid central scale ratio \f$\xi = \mu/m\f$
   * @param xi (linear) central scale ratio \f$\xi = \mu/m\f$
   */
  void setCentralScaleRatio(cdbl xi = 1.) {
    this->xi = xi;
    this->mu2 = xi * xi * this->m2;
    this->logC = 2. * log(xi);
    this->setScaleRatios(this->xiR, this->xiF);
  }

  /**
   * @brief Set reference PDF (and alpha_s)
   * @param setname PDF set name
   * @param member PDF member
   */
  void setPDF(const str& setname, cuint member) { this->resetPDF(LHAPDF::mkPDF(setname, member)); }

  /**
   * @brief Set reference PDF (and alpha_s)
   * @param setname PDF set name and member separated by /
   */
  void setPDF(const str& setname_nmem) { this->resetPDF(LHAPDF::mkPDF(setname_nmem)); }

  /**
   * @brief was a PDF set?
   * @return PDF set?
   */
  bool hasPDF() const { return static_cast<bool>(this->pdf); }

  /**
   * @brief has physical parameters?
   * @return has enough energy to do something?
   */
  bool isPhysical() const { return this->S_h > 4. * this->m2; }

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
    if (fabs(this->v.rho - 1.) < 1e-10) {
      f[0] = 0.;
      return;
    }
    // Collect all pieces
    dbl tot = 0.;
    // fill any channel
#define addLumiChannel(BIG, small)                                                        \
  if ((this->lumi_mask & LUMI_##BIG) == LUMI_##BIG) {                                     \
    tot += this->fillLumi(this->flux_##small(), this->IDX_LUMI_##BIG, FixedOrder::small); \
  }
    // gluon-gluon channel
    addLumiChannel(GG, gg);
    // quark-antiquark channel
    addLumiChannel(QQBAR, qqbar);
    // gluon-quark channel
    addLumiChannel(GQ, gq);
    // quark-quark channel
    addLumiChannel(QQ, qq);
    // quark-antiquark' channel
    addLumiChannel(QQBARPRIME, qqbarprime);
    // quark-quark' channel
    addLumiChannel(QQPRIME, qqprime);
    f[0] = tot;
    if (!std::isfinite(tot)) throw std::domain_error("Integrand is not finite!");
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
    // quark-quark channel
    {
      std::vector<PineAPPL::LumiEntry> qq;
      for (uint pid = 1; pid <= this->nl; ++pid) {
        qq.push_back({static_cast<int>(pid), static_cast<int>(pid), 1.0});
        qq.push_back({static_cast<int>(-pid), static_cast<int>(-pid), 1.0});
      }
      lumi.add(qq);
    }
    // quark-antiquark' channel
    {
      std::vector<PineAPPL::LumiEntry> qqbarprime;
      for (uint pid1 = 1; pid1 <= this->nl; ++pid1) {
        for (uint pid2 = 1; pid2 <= this->nl; ++pid2) {
          if (pid1 == pid2) continue;
          qqbarprime.push_back({static_cast<int>(pid1), static_cast<int>(-pid2), 1.0});
          qqbarprime.push_back({static_cast<int>(-pid1), static_cast<int>(pid2), 1.0});
        }
      }
      lumi.add(qqbarprime);
    }
    // quark-quark' channel
    {
      std::vector<PineAPPL::LumiEntry> qqprime;
      for (uint pid1 = 1; pid1 <= this->nl; ++pid1) {
        for (uint pid2 = 1; pid2 <= this->nl; ++pid2) {
          if (pid1 == pid2) continue;
          qqprime.push_back({static_cast<int>(pid1), static_cast<int>(pid2), 1.0});
          qqprime.push_back({static_cast<int>(-pid1), static_cast<int>(-pid2), 1.0});
        }
      }
      lumi.add(qqprime);
    }

    std::vector<PineAPPL::Order> orders = {
        // LO
        PineAPPL::Order{2, 0, 0, 0},
        // NLO
        PineAPPL::Order{3, 0, 0, 0}, PineAPPL::Order{3, 0, 1, 0}, PineAPPL::Order{3, 0, 0, 1},
        // NNLO
        PineAPPL::Order{4, 0, 0, 0}, PineAPPL::Order{4, 0, 1, 0}, PineAPPL::Order{4, 0, 0, 1},
        PineAPPL::Order{4, 0, 2, 0}, PineAPPL::Order{4, 0, 1, 1}, PineAPPL::Order{4, 0, 0, 2}};

    // fully-inclusive cross-section
    std::vector<double> bins = {0.0, 1.0};

    // create the PineAPPL grid
    PineAPPL::KeyVal kv;
    kv.set_double("x_min", this->rho_h);
    kv.set_double("q2_min", this->mu2 * 0.99);
    kv.set_double("q2_max", this->mu2 * 1.01);
    kv.set_double("q2_bins", 1);
    this->grid.reset(new PineAPPL::Grid(lumi, orders, bins, kv));
  }

  /** @brief Optimize grid */
  void optimizeGrid() const { this->grid->optimize(); }

  /**
   * @brief add scoped metadata to grid
   * @param key key
   * @param value value
   */
  void addScopedMetadata(str key, str value) const {
#define kKernelKeyStrSize 100
    char buffer[kKernelKeyStrSize];
    snprintf(buffer, kKernelKeyStrSize, "MaunaKea::%s", key.c_str());
    this->grid->set_key_value(buffer, value);
  }

  /** @brief add local metadata to grid */
  void addLocalMetadata() const {
#define kKernelValStrSize 100
    char buffer[kKernelValStrSize];
    snprintf(buffer, kKernelValStrSize, "%e", this->m2);
    this->addScopedMetadata("m2", buffer);
    snprintf(buffer, kKernelValStrSize, "%u", this->nl);
    this->addScopedMetadata("nl", buffer);
    snprintf(buffer, kKernelValStrSize, "%u", this->order_mask);
    this->addScopedMetadata("order_mask", buffer);
    snprintf(buffer, kKernelValStrSize, "%u", this->lumi_mask);
    this->addScopedMetadata("lumi_mask", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->S_h);
    this->addScopedMetadata("hadronicS", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->xi);
    this->addScopedMetadata("xi", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->xiR);
    this->addScopedMetadata("xiR", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->xiF);
    this->addScopedMetadata("xiF", buffer);
    snprintf(buffer, kKernelValStrSize, "%s/%d", this->pdf->set().name().c_str(), this->pdf->memberID());
    this->addScopedMetadata("PDF", buffer);
    this->grid->set_key_value("y_label", "sigma_tot");
    this->grid->set_key_value("y_label_tex", "\\sigma_{tot}");
    this->grid->set_key_value("y_unit", "pb");
  }

  /**
   * @brief Write grid to disk
   * @param fp file path
   */
  void writeGrid(const str& fp) const { this->grid->write(fp); }

  /** @see MC initializer */
  void Dvegas_init() const { this->grid->scale(0.); }

  /** @see MC finalizer */
  void Dvegas_final(cuint iterations) const { this->grid->scale(1. / iterations); }
};
}  // namespace MaunaKea

#endif  // MAUNAKEA_KERNEL_HPP_
