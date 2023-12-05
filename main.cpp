#include <LHAPDF/LHAPDF.h>
//#include <dvegas/dvegas.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
#include <vector>

#include "PineAPPL.hpp"
#include "Integration.hpp"

struct IntKernel : public HepSource::Integrand {
    PineAPPL::Grid* grid;
    LHAPDF::PDF* pdf;
    double rho_h;
    double m2;
    //IntKernel(): HepSource::Integrand(2,1) {}
    //virtual ~IntKernel() {}
    void operator()(const double x[], const int k[], cdbl& weight, cdbl aux[], double f[]) {
      cdbl beta_min = 0.;
      cdbl beta_h = sqrt(1. - rho_h);
      cdbl beta_max = beta_h;
      cdbl beta = beta_min + (beta_max - beta_min) * x[0];
      cdbl rho = 1 - beta*beta;
      cdbl tau = rho_h / rho;
      cdbl x1_min = tau;
      cdbl x1_max = 1.;
      cdbl x1 = x1_min + (x1_max - x1_min) * x[1];
      cdbl x2 = tau / x1;
      cdbl raw_jac = 2.* beta/rho / x1;
      cdbl mu2 = this->m2;
      dbl w = raw_jac;
      this->grid->fill(x1, x2, mu2, 0, 0.5, 0, w*weight);
      w *= this->pdf->xfxQ2(21, x1, mu2) / x1 * this->pdf->xfxQ2(21, x2, mu2) / x2;
      f[0] = w;
    }
    void Dvegas_init() const {
      this->grid->scale(0.);
    }
    void Dvegas_final(cuint iterations) const {
      this->grid->scale(1./iterations);
    };
};

// double ker (double x, 	     void * par) 
// {
//     KerArgs * args = (KerArgs*) par;
//   /*FOcomponents * p = (FOcomponents*) par;
//   double f;
  
//   if ((p->partonchannel)=="qqbar")
//     {
//       f = (p->Flux).fluxqqbar(x) * (p->FO).FOqqbar((p->rho)/x) ;
//     }
//   else if ((p->partonchannel)=="gg")
//     {
//       f = (p->Flux).fluxgg(x) * (p->FO).FOgg((p->rho)/x) ;
//     }
//   else if ((p->partonchannel)=="qg")
//     {
//       f = (p->Flux).fluxgq(x) * (p->FO).FOgq((p->rho)/x) ;
//     }
//   else if ((p->partonchannel)=="qq")
//     {
//       f = (p->Flux).fluxqq(x) * (p->FO).FOqq((p->rho)/x) ;
//     }
//   else if ((p->partonchannel)=="qqprime")
//     {
//       f = (p->Flux).fluxqqprime(x) * (p->FO).FOqqprime((p->rho)/x) ;
//     }
//   else if ((p->partonchannel)=="qqbarprime")
//     {
//       f = (p->Flux).fluxqqbarprime(x) * (p->FO).FOqqbarprime((p->rho)/x) ;
//     }
//   else//all channels combined
//     {
//       f = (p->Flux).fluxqqbar(x)    * (p->FO).FOqqbar((p->rho)/x)      +
// 	(p->Flux).fluxgg(x)         * (p->FO).FOgg((p->rho)/x)         +
// 	(p->Flux).fluxgq(x)         * (p->FO).FOgq((p->rho)/x)         +
// 	(p->Flux).fluxqq(x)         * (p->FO).FOqq((p->rho)/x)         +
// 	(p->Flux).fluxqqprime(x)    * (p->FO).FOqqprime((p->rho)/x)    +
// 	(p->Flux).fluxqqbarprime(x) * (p->FO).FOqqbarprime((p->rho)/x) ;
//     }
//   return f/x;*/
// }

void fill_grid(PineAPPL::Grid* grid, LHAPDF::PDF* pdf, cdbl rho_h, cdbl m2) {
    IntKernel k;
    k.grid = grid;
    k.pdf = pdf;
    k.rho_h = rho_h;
    k.m2 = m2;
    IntegrationConfig cfg;
    cfg.MC_warmupCalls = 1000;
    cfg.calls = 20000;
    cfg.verbosity = 3;
    IntegrationOutput out;
    integrate2D(&k, cfg, &out);
    printf("res: %e +- %e\n", out.result, out.error);
}

int main() {
    cdbl S_h = pow(8e3,2);
    cdbl m2 = pow(172.5,2);
    cdbl rho_h = 4.*m2/S_h;
    // create a new luminosity function
    PineAPPL::Lumi lumi;
    lumi.add({PineAPPL::LumiEntry {21,21,1.0}});

    std::unique_ptr<LHAPDF::PDF> pdf (LHAPDF::mkPDF("gonly", 0));

    // only LO
    std::vector<PineAPPL::Order> orders = {PineAPPL::Order {0,0,0,0}};

    // fully-inclusive cross-section
    std::vector<double> bins = {0.0, 1.0};

    // create the PineAPPL grid with default interpolation and binning parameters
    PineAPPL::KeyVal kv;
    PineAPPL::Grid grid(lumi, orders, bins, kv);

    // fill the grid with phase-space points
    fill_grid(&grid, pdf.get(), rho_h, m2);
    grid.optimize();

    // store some metadata in the grid
    //grid.set_key_value("events", "10000000");

    // write the grid to disk - with `.lz4` suffix the grid is automatically LZ4 compressed
    const std::string filename = "MaunaKea.pineappl.lz4";
    grid.write(filename);
}
