#ifndef FO_HPP_
#define FO_HPP_

#include <gsl/gsl_sf.h>

#include "./config.h"

// some global constants
#define pi M_PI
#define CA 3.0
#define TR 0.5

namespace MaunaKea {

/** @name LO results */
///@{
dbl f0gg(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl beta = sqrt(1 - rho);
  cdbl f =
      pi * beta * (rho / 192.) * ((pow(rho, 2) + 16 * rho + 16) * log((1 + beta) / (1 - beta)) / beta - 28 - 31 * rho);
  return f;
}

dbl f0qqbar(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl beta = sqrt(1 - rho);
  return pi * beta * rho / 27. * (2 + rho);
}
///@}

/** @name NLO results */
///@{

//==========================================================================
//   High-quality parameterization exact NLO results (Czakon, Mitov '08)
//   Fits from HATHOR: Aliev et al arXiv:1007.1327
//==========================================================================

dbl f1qqbar(cdbl rho, cdbl nl) {
  // double nl = 5.0;
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lbe2 = pow(lbe, 2);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  cdbl rho2 = pow(rho, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);

  cdbl f10nl0 = +0.894801363487060697 * b2 - 15.9806533566333508 * b3 + 15.5948076822314514 * b4 -
                0.508993720940668306 * b5 + 0.25808023534387927 * b2 * lbe - 3.99149865104188887 * b3 * lbe -
                8.97737571174816865 * b4 * lbe + 0.147026822789739241 * b5 * lbe + 0.0187135812730567129 * b2 * lbe2 -
                1.81602862945435741 * b3 * lbe2 - 1.74752529849853127 * b * lrho + 6.20751624259016067 * b2 * lrho -
                6.75340649699489205 * b3 * lrho + 2.29340856664909682 * b4 * lrho + 0.135309206882346187 * b * lrho2 -
                0.0712992055892768895 * b2 * lrho2 - 0.0640103309259597209 * b3 * lrho2 - 0.0913852259360125798 * rho +
                0.44171954644071765 * b * rho - 0.57251372837945371 * b * lbe * rho +
                1.185185185185185185 * b * lbe2 * rho + 0.0987654320987654321 * b * lrho * rho +
                0.0138455459529272122 * b * rho2 + 0.049382716049382716 * b * lrho * rho2;

  cdbl f10nl1 = (b * (3 - 4 * b2 + b4) * (-5 - 3 * lrho + log(64))) / 243.;
  return (f10nl0 + nl * f10nl1);
}

dbl f1gg(cdbl rho, cdbl nl) {
  // double nl = 5.0;
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lbe2 = pow(lbe, 2);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  cdbl rho2 = pow(rho, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);
  cdbl logxconf = log((1.0 - b) / (1.0 + b));

  cdbl f10nl0 = +0.0471205071232564865 + 0.235655553965240606 * b - 112.1628025020790127 * b2 +
                1883.77093917245376 * b3 - 1766.27896687670163 * b4 - 4.28709466299521533 * b5 -
                0.086738651030143563 * b * lbe - 30.2912153871041909 * b2 * lbe + 687.805453265266195 * b3 * lbe +
                1142.477618105764338 * b4 * lbe - 61.3742807324254932 * b5 * lbe + 0.875 * b * lbe2 -
                2.19494021495237345 * b2 * lbe2 + 169.273238020325875 * b3 * lbe2 + 284.814617211531812 * b * lrho -
                849.780999197478008 * b2 * lrho + 817.825794348102379 * b3 * lrho - 252.863693983687393 * b4 * lrho +
                57.8966224150269576 * b * lrho2 - 121.9429831242445605 * b2 * lrho2 + 64.0461741451940366 * b3 * lrho2 +
                0.03125 * b * rho2 + 0.015625 * logxconf * rho2;

  cdbl f10nl1 = -((2 * b + logxconf) * rho2) / 256.;
  return (f10nl0 + nl * f10nl1);
}

dbl f1gq(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  // cdbl lbe2 = pow(lbe, 2);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  cdbl lrho3 = lrho * lrho2;
  cdbl lrho4 = pow(lrho2, 2);
  cdbl lrho5 = lrho * lrho4;
  // cdbl rho2 = pow(rho, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);

  cdbl f10nl0 = -0.02457581886482620869 * b3 - 3.28032159897403943 * b4 + 3.7941231153033735 * b5 -
                0.189185108551928629 * b6 + 0.1388888888888888889 * b3 * lbe - 0.0178630821547912901 * b4 * lbe -
                0.585680726926105217 * b5 * lbe - 1.8961443995961806 * b6 * lbe -
                3.19157673250705547 * b2 * lrho * rho + 5.01130734325825868 * b3 * lrho * rho -
                1.78170350247765988 * b4 * lrho * rho - 0.1255542812553638502 * b2 * lrho2 * rho -
                0.307143757144090147 * b3 * lrho2 * rho + 0.23465015682704749 * b4 * lrho2 * rho +
                0.0299903911910146112 * b2 * lrho3 * rho - 0.000427110882824291123 * b2 * lrho4 * rho -
                0.00001115993665476662179 * b2 * lrho5 * rho;
  return f10nl0;
}

//==========================================================================
//   Scale-dependent terms (as in: Nason, Dawson, Ellis '88)
//==========================================================================

dbl beta0(cdbl nl) { return (11. / 3. * CA - 4. / 3. * TR * nl) / (4. * pi); }

dbl fbarR1qqbar(cdbl rho, cdbl nl) {
  cdbl f = 2. * beta0(nl) * f0qqbar(rho, nl);
  return f;
}

dbl fbarR1gg(cdbl rho, cdbl nl) {
  cdbl f = 2. * beta0(nl) * f0gg(rho, nl);
  return f;
}

dbl h1(cdbl beta) {
  cdbl f = -2 * gsl_sf_dilog((1 - beta) / 2.) + 2 * gsl_sf_dilog((1 + beta) / 2.) - pow(log((1 - beta) / 2.), 2) +
           pow(log((1 + beta) / 2.), 2);
  return f;
}

dbl h2(cdbl beta) {
  cdbl f = -gsl_sf_dilog((-2 * beta) / (1 - beta)) + gsl_sf_dilog((2 * beta) / (1 + beta));
  return f;
}

dbl fbar1qqbar(cdbl rho, cdbl nl) {
  cdbl beta = sqrt(1 - rho);
  // double nl = 5;
  cdbl f = (2 * rho * log((1 + beta) / (1 - beta))) / (81. * pi) +
           f0qqbar(rho, nl) *
               (-(-127 + 6 * nl) / (72. * pow(pi, 2)) + (2 * log(rho / (4. * pow(beta, 2)))) / (3. * pow(pi, 2)));
  return f * 4. * pi;
}

dbl fbarF1qqbar(cdbl rho, cdbl nl) { return fbar1qqbar(rho, nl) - fbarR1qqbar(rho, nl); }

dbl fbar1gg(cdbl rho, cdbl nl) {
  cdbl beta = sqrt(1 - rho);
  cdbl f = (-181 * beta) / (1440. * pi) + (26 * beta * rho) / (45. * pi) - (2483 * beta * pow(rho, 2)) / (1920. * pi) +
           (-rho / (8. * pi) + pow(rho, 2) / (16. * pi) - pow(rho, 3) / (256. * pi)) * h1(beta) +
           (rho / (8. * pi) + pow(rho, 2) / (8. * pi) + pow(rho, 3) / (128. * pi)) * h2(beta) +
           ((-3 * rho) / (8. * pi) + (33 * pow(rho, 2)) / (128. * pi) + (59 * pow(rho, 3)) / (768. * pi)) *
               log((1 + beta) / (1 - beta)) +
           (3 * f0gg(rho, nl) * log(rho / (4. * pow(beta, 2)))) / (2. * pow(pi, 2));
  return f * 4. * pi;
}

dbl fbarF1gg(cdbl rho, cdbl nl) { return fbar1gg(rho, nl) - fbarR1gg(rho, nl); }

dbl fbarF1gq(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl beta = sqrt(1 - rho);
  cdbl f = (-181 * beta) / (6480. * pi) + (289 * beta * rho) / (2160. * pi) -
           (1319 * beta * pow(rho, 2)) / (25920. * pi) + (-rho / (72. * pi) + pow(rho, 2) / (144. * pi)) * h1(beta) +
           ((-17 * rho) / (432. * pi) + pow(rho, 2) / (128. * pi) + (7 * pow(rho, 3)) / (1728. * pi)) *
               log((1 + beta) / (1 - beta));
  return f * 4. * pi;
}

///@}

/** @name Coefficient function collections */
///@{

/** @brief Single coefficient function */
typedef dbl (*coeff)(cdbl rho, cdbl nl);

/** @brief Map of coefficient functions */
struct CoeffMap {
  /** @brief LO */
  coeff f0;
  /** @brief NLO */
  coeff f1;
  /** @brief NLO fact. SV */
  coeff fbarF1;
  /** @brief NLO ren. SV */
  coeff fbarR1;
};

/** @brief gluon-gluon channel */
const CoeffMap gg = {f0gg, f1gg, fbarF1gg, fbarR1gg};

/** @brief quark-antiquark channel */
const CoeffMap qqbar = {f0qqbar, f1qqbar, fbarF1qqbar, fbarR1qqbar};

/** @brief gluon-quark channel */
const CoeffMap gq = {0, f1gq, fbarF1gq, 0};
///@}

}  // namespace MaunaKea

#endif  // FO_HPP_
