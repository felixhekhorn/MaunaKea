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
  /** @brief renormalization scale ratio \f$\xi_R = \mu_R/m\f$ */
  dbl xiR;
  /** @brief renormalization scale log \f$L_R = 2*\log(\xi_R)\f$ */
  dbl logR;
  /** @brief strong coupling constant */
  dbl as;
  /** @brief factorization scale ratio \f$\xi_F = \mu_F/m\f$ */
  dbl xiF;
  /** @brief factorization scale log \f$L_F = 2*\log(\xi_F)\f$ */
  dbl logF;
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
   * @brief Insert a given luminosity into the grid
   * @param flux flux
   * @param idx_lumi grid luminosity index
   * @param m coefficient functions
   */
  dbl fillLumi(cdbl flux, cuint idx_lumi, const FixedOrder::CoeffMap m) const {
    dbl tot = 0.;
    // fill any order
#define fillOrder(label, k, idx, powR, powF)                                                     \
  if (m.f##label##k) {                                                                           \
    cdbl weight = this->v.common_weight * m.f##label##k(this->v.rho, this->nl);                  \
    this->grid->fill(this->v.x1, this->v.x2, this->m2, this->IDX_ORDER_##idx, 0.5, idx_lumi,     \
                     weight* this->v.vegas_weight* this->v.x1* this->v.x2);                      \
    tot += weight * flux * pow(this->as, 2 + k) * pow(this->logR, powR) * pow(this->logF, powF); \
  }
    // LO
    if ((this->order_mask & ORDER_LO) == ORDER_LO) {
      fillOrder(, 0, LO, 0, 0);
    }
    // NLO
    if ((this->order_mask & ORDER_NLO) == ORDER_NLO) {
      fillOrder(, 1, NLO, 0, 0);
      // SV
      fillOrder(barR, 1, NLO_R, 1, 0);
      fillOrder(barF, 1, NLO_F, 0, 1);
    }
    // NNLO
    if ((this->order_mask & ORDER_NNLO) == ORDER_NNLO) {
      fillOrder(, 2, NNLO, 0, 0);
      // SV^1
      fillOrder(barR, 2, NNLO_R, 1, 0);
      fillOrder(barF, 2, NNLO_F, 0, 1);
      // SV^2
      fillOrder(barRR, 2, NNLO_RR, 2, 0);
      fillOrder(barRF, 2, NNLO_RF, 1, 1);
      fillOrder(barFF, 2, NNLO_FF, 0, 2);
    }
    return tot;
  }

  /** @brief Cache alpha_s at the renormalization scale */
  void setAlphaS() {
    cdbl muR2 = this->xiR * this->xiR * this->m2;
    this->as = this->pdf->alphasQ2(muR2);
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
  /** @brief quark-quark marker */
  static cuint LUMI_QQ = 1 << 3;
  /** @brief quark-antiquark' marker */
  static cuint LUMI_QQBARPRIME = 1 << 4;
  /** @brief quark-quark' marker */
  static cuint LUMI_QQPRIME = 1 << 5;
  /** @brief all order marker */
  static cuint LUMI_ALL = LUMI_GG | LUMI_QQBAR | LUMI_GQ | LUMI_QQ | LUMI_QQBARPRIME | LUMI_QQPRIME;
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
   * @brief Set renormalization and factorization scale ratio \f$\xi = \mu/m\f$
   * @param xiR renormalization scale ratio \f$\xi_R = \mu_R/m\f$
   * @param xiF factorization scale ratio \f$\xi_R = \mu_F/m\f$
   */
  void setScaleRatios(cdbl xiR = 1., cdbl xiF = 1.) {
    this->xiR = xiR;
    this->logR = 2. * log(xiR);
    // recalculate as
    if (this->hasPDF()) this->setAlphaS();
    this->xiF = xiF;
    this->muF2 = xiF * xiF * this->m2;
    this->logF = 2. * log(xiF);
  }

  /**
   * @brief Set reference PDF (and alpha_s)
   * @param setname PDF set name
   * @param member PDF member
   */
  void setPDF(const str setname, cuint member) {
    this->pdf.reset(LHAPDF::mkPDF(setname, member));
    this->setAlphaS();
  }

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
    // Collect all pieces
    dbl tot = 0.;
    // fill any channel
#define addLumiChannel(BIG, small)                                                        \
  if ((this->lumi_mask & this->LUMI_##BIG) == this->LUMI_##BIG) {                         \
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
    kv.set_double("q2_min", this->m2 * 0.99);
    kv.set_double("q2_max", this->m2 * 1.01);
    kv.set_double("q2_bins", 1);
    this->grid.reset(new PineAPPL::Grid(lumi, orders, bins, kv));
  }

  /** @brief Optimize grid */
  void optimizeGrid() const { this->grid->optimize(); }

  /**
   * @brief add raw metadata to grid
   * @param key key
   * @param value value
   */
  void addRawMetadata(str key, str value) const {
#define kKernelKeyStrSize 100
    char buffer[kKernelKeyStrSize];
    snprintf(buffer, kKernelKeyStrSize, "MaunaKea/%s", key.c_str());
    this->grid->set_key_value(buffer, value);
  }

  /**
   * @brief add local metadata to grid
   */
  void addLocalMetadata() const {
#define kKernelValStrSize 100
    char buffer[kKernelValStrSize];
    snprintf(buffer, kKernelValStrSize, "%e", this->m2);
    this->addRawMetadata("m2", buffer);
    snprintf(buffer, kKernelValStrSize, "%u", this->nl);
    this->addRawMetadata("nl", buffer);
    snprintf(buffer, kKernelValStrSize, "%u", this->order_mask);
    this->addRawMetadata("order_mask", buffer);
    snprintf(buffer, kKernelValStrSize, "%u", this->lumi_mask);
    this->addRawMetadata("lumi_mask", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->S_h);
    this->addRawMetadata("hadronicS", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->xiR);
    this->addRawMetadata("xiR", buffer);
    snprintf(buffer, kKernelValStrSize, "%e", this->xiF);
    this->addRawMetadata("xiF", buffer);
    snprintf(buffer, kKernelValStrSize, "%s#%d", this->pdf->set().name().c_str(), this->pdf->memberID());
    this->addRawMetadata("PDF", buffer);
    this->grid->set_key_value("y_label", "sigma_tot");
    this->grid->set_key_value("y_label_tex", "\\sigma_{tot}");
    this->grid->set_key_value("y_label", "pb");
  }

  /**
   * @brief Write grid to disk
   * @param fp file path
   */
  void writeGrid(const str fp) const { this->grid->write(fp); }

  /** @see MC initializer */
  void Dvegas_init() const { this->grid->scale(0.); }

  /** @see MC finalizer */
  void Dvegas_final(cuint iterations) const { this->grid->scale(1. / iterations); }
};
}  // namespace MaunaKea

#endif  // KERNEL_HPP_
