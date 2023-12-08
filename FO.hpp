#ifndef FO_HPP_
#define FO_HPP_

#include "./config.h"

namespace MaunaKea {

/** @name LO results */
///@{
dbl f0gg(cdbl rho) {
  cdbl beta = sqrt(1 - rho);
  cdbl f = M_PI * beta * (rho / 192.) *
           ((pow(rho, 2) + 16 * rho + 16) * log((1 + beta) / (1 - beta)) / beta - 28 - 31 * rho);
  return f;
}

dbl f0qqbar(cdbl rho) {
  cdbl beta = sqrt(1 - rho);
  return M_PI * beta * rho / 27. * (2 + rho);
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

dbl f1gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lbe2 = pow(lbe, 2);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  cdbl lrho3 = lrho * lrho2;
  cdbl lrho4 = pow(lrho2, 2);
  cdbl lrho5 = lrho * lrho4;
  cdbl rho2 = pow(rho, 2);
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

///@}

}  // namespace MaunaKea

#endif  // FO_HPP_
