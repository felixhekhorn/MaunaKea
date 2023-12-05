#include <LHAPDF/LHAPDF.h>
//#include <dvegas/dvegas.h>

#include <cmath>
#include <cstddef>
#include <cstdio>
#include <random>
#include <vector>

#include "PineAPPL.hpp"
//#include "Integration.hpp"

struct IntKernel /*: public HepSource::Integrand*/ {
    PineAPPL::Grid* grid;
    LHAPDF::PDF* pdf;
    double rho_h;
    double m2;
    //IntKernel(): HepSource::Integrand(2,1) {}
    //virtual ~IntKernel() {}
    void operator()(const double x[], const int k[], const double& weight, const double aux[], double f[]) {
      const double beta_min = 0.;
      const double beta_h = sqrt(1. - rho_h);
      const double beta_max = beta_h;
      const double beta = beta_min + (beta_max - beta_min) * x[0];
      const double rho = 1 - beta*beta;
      const double tau = rho_h / rho;
      const double x1_min = tau;
      const double x1_max = 1.;
      // const double x1 = x1_min + (x1_max - x1_min) * x[1];
      // const double x2 = tau / x1;
      double x1 = x[0];
      double x2 = x[1];
      const double raw_jac = 2.* beta/rho / x1;
      const double mu2 = this->m2;
      double w = raw_jac;
      w = 1.;
      w *= weight;
      this->grid->fill(x1, x2, mu2, 0, 0.5, 0, w);
      w *= this->pdf->xfxQ2(21, x1, mu2) * this->pdf->xfxQ2(21, x2, mu2);
      f[0] = w;
    }
    //void Dvegas_init() const {}
    //void Dvegas_final(cuint iterations) const {};
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

double fill_grid(PineAPPL::Grid* grid, LHAPDF::PDF* pdf, double rho_h) {
    // IntKernel k;
    // k.grid = &grid;
    // k.pdf = &pdf;
    // k.rho_h = rho_h;
    // k.m2 = pow(172.5,2);
}


// double fill_grid(PineAPPL::Grid &grid, LHAPDF::PDF &pdf, double rho_h) {
//     IntKernel k;
//     k.grid = &grid;
//     k.pdf = &pdf;
//     k.rho_h = rho_h;
//     k.m2 = pow(172.5,2);

//     double res, err;

//     double xl[2] = { 0., 0.};
//     double xu[2] = { 1., 1.};

//     const gsl_rng_type *T;
//     gsl_rng *r;

//     gsl_monte_function G = { &callFunctor2D<IntKernel>, 2, &k };

//     size_t calls = 100000;

//     gsl_rng_env_setup ();

//     T = gsl_rng_default;
//     r = gsl_rng_alloc (T);
//     gsl_monte_vegas_state *s = gsl_monte_vegas_alloc (2);

//     gsl_monte_vegas_integrate (&G, xl, xu, 2, 10000, r, s,
//                                &res, &err);
//     printf ("vegas warm-up %e +- %e\n", res, err);

//     printf ("converging...\n");

//     unsigned int guard = 0;
//     do
//       {
//         gsl_monte_vegas_integrate (&G, xl, xu, 2, calls, r, s, &res, &err);
//         printf ("result = %e sigma = %e "
//                 "chisq/dof = %.1f\n", res, err, gsl_monte_vegas_chisq (s));
//       }
//     while (fabs (gsl_monte_vegas_chisq (s) - 1.0) > 0.5 && guard++ < 20);

//     printf ("vegas final %e +- %e\n", res, err);

//     gsl_monte_vegas_free (s);

//   gsl_rng_free (r);
//   return res;
// }

int main() {
    // create a new luminosity function for the $\gamma\gamma$ initial state
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
    double res = fill_grid(&grid, pdf.get(), 0.1);

    // store some metadata in the grid
    //grid.set_key_value("events", "10000000");

    // write the grid to disk - with `.lz4` suffix the grid is automatically LZ4 compressed
    const std::string filename = "MaunaKea.pineappl.lz4";
    grid.write(filename);
}
