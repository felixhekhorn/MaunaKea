#ifndef FO_HPP_
#define FO_HPP_

#include <gsl/gsl_sf.h>

#include "./config.h"

// some global constants
#define pi M_PI
#define CA 3.0
#define TR 0.5

namespace MaunaKea {

/** @brief Fixed order coefficient functions */
namespace FixedOrder {

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
//   Scale-dependent terms (as in: Nason, Dawson, Ellis '88) (recall to rescale!)
//==========================================================================

dbl beta0(cdbl nl) { return (11. / 3. * CA - 4. / 3. * TR * nl) / (4. * pi); }

dbl fbarR1qqbar(cdbl rho, cdbl nl) {
  cdbl f = 2. * beta0(nl) * f0qqbar(rho, dblNaN);
  return f;
}

dbl fbarR1gg(cdbl rho, cdbl nl) {
  cdbl f = 2. * beta0(nl) * f0gg(rho, dblNaN);
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
           f0qqbar(rho, dblNaN) *
               (-(-127 + 6 * nl) / (72. * pow(pi, 2)) + (2 * log(rho / (4. * pow(beta, 2)))) / (3. * pow(pi, 2)));
  return f * 4. * pi;
}

dbl fbarF1qqbar(cdbl rho, cdbl nl) { return fbar1qqbar(rho, nl) - fbarR1qqbar(rho, nl); }

dbl fbar1gg(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl beta = sqrt(1 - rho);
  cdbl f = (-181 * beta) / (1440. * pi) + (26 * beta * rho) / (45. * pi) - (2483 * beta * pow(rho, 2)) / (1920. * pi) +
           (-rho / (8. * pi) + pow(rho, 2) / (16. * pi) - pow(rho, 3) / (256. * pi)) * h1(beta) +
           (rho / (8. * pi) + pow(rho, 2) / (8. * pi) + pow(rho, 3) / (128. * pi)) * h2(beta) +
           ((-3 * rho) / (8. * pi) + (33 * pow(rho, 2)) / (128. * pi) + (59 * pow(rho, 3)) / (768. * pi)) *
               log((1 + beta) / (1 - beta)) +
           (3 * f0gg(rho, dblNaN) * log(rho / (4. * pow(beta, 2)))) / (2. * pow(pi, 2));
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

/** @name NNLO results */
///@{

//=========================================================================
//   threshold expansion of the NNLO function f20
//   (with exact BORN factored out) from Beneke et al (2009); C2gg=0.
//   The constant C2gg is implicitly contained in the fits below
//   (but note it also receives contrubution from: Born x 1/beta^2)
//=========================================================================

dbl PXS_NNLO_THRESHOLD_gg(cdbl rho, cdbl nl) {
  // double nl = 5.0;// number of light flavors
  cdbl b = sqrt(1 - rho);
  cdbl b2 = pow(b, 2);
  cdbl lbe = log(b);
  cdbl lbe2 = pow(lbe, 2);
  cdbl lbe3 = pow(lbe, 3);
  cdbl lbe4 = pow(lbe, 4);

  cdbl thresh20nl0 = -0.024007157914810687 / b + 0.43407982319605976 / b2 + 14.861847493860669 * lbe +
                     (1.8153673293788755 * lbe) / b - 1.9983842014515505 * lbe2 + (3.142857142857143 * lbe2) / b -
                     14.701583665657568 * lbe3 + 29.18050088899328 * lbe4;

  cdbl thresh20nl1 = -0.0061192368274098 / b + 0.13912366816397184 * lbe + (0.04365079365079365 * lbe) / b -
                     0.7558262174090832 * lbe2 + 0.5403796460924681 * lbe3;

  cdbl f = thresh20nl0 + nl * thresh20nl1;
  return f0gg(rho, dblNaN) * f;
}

//===========================================================================
//    The power suppressed terms (no scales).
//    Includes the gg "constant" term.
//===========================================================================

dbl f20nl2_FIT_gg(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  // cdbl lbe = log(b);
  cdbl rho2 = pow(rho, 2);
  cdbl rho3 = pow(rho, 3);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);

  cdbl f = (6.44022 * b * rho - 4.8664 * b2 * rho - 0.0324653 * lrho2 * rho - 13.8424 * b * rho2 + 4.7366 * b2 * rho2 -
            2.91398 * lrho * rho2 + 8.43828 * b * rho3 - 2.78748 * b2 * rho3 + 2.38971 * b3 * rho3) /
           10000.0;
  return f;
}

dbl f20nl1_FIT_gg(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lbe2 = pow(lbe, 2);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  cdbl lrho3 = pow(lrho, 3);
  cdbl rho2 = pow(rho, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);
  cdbl b7 = pow(b, 7);

  cdbl f = -0.0195046 * b - 1.4717 * b2 - 0.223616 * b3 + 0.499196 * b5 + 1.32756 * b7 + 0.00466872 * b3 * lbe +
           0.0321469 * b6 * lbe2 + 0.579781 * lrho2 * rho + 0.166646 * lrho3 * rho - 1.36644 * lrho * rho2 +
           2.24909 * lrho2 * rho2;

  return f;
}

dbl f20nl0_FIT_gg(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lrho = log(rho);
  cdbl lrho2 = pow(lrho, 2);
  cdbl rho2 = pow(rho, 2);
  cdbl rho3 = pow(rho, 3);
  cdbl rho4 = pow(rho, 4);
  cdbl rho5 = pow(rho, 5);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);

  cdbl f = 581.27542 * b + 1251.4057 * b2 - 60.478096 * b3 + 1101.2272 * b4 - 2905.3858 * b5 + 629.9128 * b4 * lbe -
           5.189107503169016 * lrho + 1200.741 * lrho * rho + 162.50333 * lrho2 * rho - 1810.2849 * b * rho2 +
           36.074524 * lrho * rho2 - 1192.8918 * lrho2 * rho2 + 1568.7591 * b * rho3 - 461.21326 * b * rho4 +
           121.6379 * b * rho5;

  return f;
}

dbl PXS_NNLO_FIT_gg(cdbl rho, cdbl nl) {
  // double nl=5.0;
  cdbl nl2 = pow(nl, 2);
  // the functions fij are normalized in as^n/m^2
  cdbl f20 = nl2 * f20nl2_FIT_gg(rho) + nl * f20nl1_FIT_gg(rho) + f20nl0_FIT_gg(rho);
  return f20;
}

dbl f2gg(cdbl rho, cdbl nl) {
  cdbl f = PXS_NNLO_THRESHOLD_gg(rho, nl) + PXS_NNLO_FIT_gg(rho, nl);
  return f;
}

//================================================================================
// The part of sigma(q qbar -> t tbar q qbar) not included in PXS_NNLO_FIT_qqbar.
// It is a sum of a pure interference term and sigma(q qbar' -> t tbar q qbar').
//================================================================================

dbl f20_FIT_qqbarprime(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -0.740572 - 31.2117 * pow(b, 2) - 0.31495 * pow(b, 3) + 15.8601 * pow(b, 4) - 1.64639 * pow(b, 5) +
              18.9767 * pow(b, 6) - 19.6977 * lrho * rho - 3.16565 * pow(lrho, 2) * rho - 16.1386 * lrho * pow(rho, 2) +
              12.3828 * pow(lrho, 2) * pow(rho, 2) + 4.17707 * lrho * pow(rho, 3);
  return -0.4768323995789214 * lrho - pow(b, 2) * exp(fexp);
}

dbl PXS_NNLO_identical_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lrho = log(rho);
  // the fit of the pure interference contribution:
  cdbl f = -1.26644 * pow(b, 20) * lbe + 1.53647 * pow(b, 3) * rho + 10.7411 * pow(b, 4) * rho + 3.86236 * lrho * rho +
           0.327143 * pow(lrho, 2) * rho - 24.3298 * pow(b, 4) * pow(rho, 2) - 21.332 * lrho * pow(rho, 2) -
           10.7669 * pow(lrho, 2) * pow(rho, 2) - 4.50719 * pow(b, 3) * pow(rho, 3) +
           15.4975 * pow(b, 4) * pow(rho, 3) + 17.4705 * lrho * pow(rho, 3) + 2.90068 * pow(b, 3) * pow(rho, 4) -
           4.98808 * pow(b, 4) * pow(rho, 4);

  return f + f20_FIT_qqbarprime(rho);
}

//=========================================================================
//   threshold expansion of the NNLO functions f22, f21 and f20
//   (with exact BORN factored out) from Beneke et al (2009); C2qq=0.
//   The actual constant is set to zero here and is contained in the fits.
//=========================================================================

dbl PXS_NNLO_THRESHOLD_qqbar(cdbl rho, cdbl nl) {
  // double nl = 5.0;// number of light flavors
  cdbl b = sqrt(1 - rho);
  cdbl lbe = log(b);
  cdbl lbe2 = pow(lbe, 2);
  cdbl lbe3 = pow(lbe, 3);
  cdbl lbe4 = pow(lbe, 4);

  cdbl thresh20nl0 = 0.022846306484003143 / pow(b, 2) - 0.033390549389716966 / b + 1.5810883840639216 * lbe +
                     (0.3422025061689375 * lbe) / b + 6.62693178820539 * lbe2 - (0.8888888888888888 * lbe2) / b -
                     9.531529857990602 * lbe3 + 5.76404955831966 * lbe4;
  cdbl thresh20nl1 = 0.01168217939778234 / b + 0.35320738606771657 * lbe - (0.027777777777777776 * lbe) / b -
                     0.5752392118254954 * lbe2 + 0.2401687315966525 * lbe3;
  // cdbl thresh21nl0 = 2.963036815826603 - 0.5208333333333334/b - 5.322595839415953*lbe
  //   + (0.4444444444444444*lbe)/b + 11.307833327841358*lbe2 - 5.76404955831966*lbe3;
  // cdbl thresh21nl1 = -0.40917468164069626 + 0.041666666666666664/b
  //   + 0.46164291174543004*lbe - 0.5403796460924681*lbe2;
  // cdbl thresh22nl0 = 1.7154768057697534 - 3.5187082038820767*lbe + 1.441012389579915*lbe2;
  // cdbl thresh22nl1 = -0.26328935946926935 + 0.22515818587186173*lbe;

  cdbl f20 = thresh20nl0 + nl * thresh20nl1;
  // cdbl f21 = thresh21nl0 + nl*thresh21nl1;
  // cdbl f22 = thresh22nl0 + nl*thresh22nl1;

  // cdbl f = f20 + f21*(lg + lgfr) + f22*(pow(lg,2) + 2.*lg*lgfr + pow(lgfr,2));
  return f0qqbar(rho, dblNaN) * f20;
}

//==========================================================================
//   The fits for qqbar -> ttbar:
//==========================================================================

dbl f22nl2_FIT_qqbar(cdbl rho) { return f0qqbar(rho, dblNaN) / (12. * pow(pi, 2)); }

dbl f22nl1_FIT_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl f = 0.00239658 * b * rho - 0.0220604 * pow(b, 2) * rho + 0.0188993 * pow(b, 3) * rho - 0.0131116 * lrho * rho -
           0.000963099 * b * pow(rho, 2) + 0.0141569 * pow(b, 2) * pow(rho, 2) - 0.00104753 * pow(b, 3) * pow(rho, 2) -
           0.000699993 * lrho * pow(rho, 2) - 0.0014341 * b * pow(rho, 3) + 0.00610629 * lrho * pow(rho, 3) -
           0.000214179 * lrho * pow(rho, 4);
  return f;
}

dbl f22nl0_FIT_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl f = 0.191165 * pow(b, 2) - 13.2537 * b * rho + 195.393 * pow(b, 2) * rho - 139.565 * pow(b, 3) * rho +
           10.4168 * lrho * rho + 16.2786 * b * pow(rho, 2) + 165.851 * pow(b, 2) * pow(rho, 2) +
           235.88 * pow(b, 3) * pow(rho, 2) + 213.556 * lrho * pow(rho, 2) - 3.03232 * b * pow(rho, 3) +
           835.447 * pow(b, 2) * pow(rho, 3) - 90.1477 * pow(b, 3) * pow(rho, 3) + 677.134 * lrho * pow(rho, 3) +
           295.537 * lrho * pow(rho, 4);
  return f;
}

dbl f21nl2_FIT_qqbar(cdbl rho) {
  cdbl lrho = log(rho);
  return f0qqbar(rho, dblNaN) * (5 + 3 * lrho - 6 * log(2)) / (18. * pow(pi, 2));
}

dbl f21nl1_FIT_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl f = -0.000188644 * pow(b, 2) + 0.371362 * b * rho - 3.14887 * pow(b, 2) * rho + 3.50366 * pow(b, 3) * rho +
           0.102822 * lrho * rho - 0.328962 * b * pow(rho, 2) + 8.58138 * pow(b, 2) * pow(rho, 2) -
           3.93341 * pow(b, 3) * pow(rho, 2) + 1.55362 * lrho * pow(rho, 2) - 0.0428876 * b * pow(rho, 3) -
           1.31531 * pow(b, 2) * pow(rho, 3) + 0.896763 * pow(b, 3) * pow(rho, 3) + 2.39531 * lrho * pow(rho, 3) +
           0.0208085 * lrho * pow(rho, 4);
  return f;
}

dbl f21nl0_FIT_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl f = -0.720681 * pow(b, 2) + 46.5167 * b * rho - 761.601 * pow(b, 2) * rho + 503.018 * pow(b, 3) * rho -
           50.2996 * lrho * rho - 61.2529 * b * pow(rho, 2) - 1148.19 * pow(b, 2) * pow(rho, 2) -
           914.708 * pow(b, 3) * pow(rho, 2) - 1046.69 * lrho * pow(rho, 2) + 14.8024 * b * pow(rho, 3) -
           3847.6 * pow(b, 2) * pow(rho, 3) + 377.917 * pow(b, 3) * pow(rho, 3) - 3269.96 * lrho * pow(rho, 3) -
           1386.84 * lrho * pow(rho, 4);
  return f;
}

dbl f20nl2_FIT_qqbar(cdbl rho) {
  cdbl lrho = log(rho);
  cdbl f = (25 - 3 * pow(pi, 2) + 30 * (lrho - 2 * log(2)) + 9 * pow(lrho - 2 * log(2), 2)) / (108. * pow(pi, 2));
  return f * f0qqbar(rho, dblNaN);
}

dbl f20nl1_FIT_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl f = 0.90756 * b * rho - 6.75556 * pow(b, 2) * rho + 9.18183 * pow(b, 3) * rho + 1.3894 * lrho * rho +
           0.137881 * pow(lrho, 2) * rho - 0.99749 * b * pow(rho, 2) + 27.7454 * pow(b, 2) * pow(rho, 2) -
           12.9055 * pow(b, 3) * pow(rho, 2) + 6.13693 * lrho * pow(rho, 2) - 0.0077383 * b * pow(rho, 3) -
           4.49375 * pow(b, 2) * pow(rho, 3) + 3.86854 * pow(b, 3) * pow(rho, 3) + 8.78276 * lrho * pow(rho, 3) -
           0.380386 * pow(b, 4) * pow(rho, 4) - 0.0504095 * lrho * pow(rho, 4);
  return f;
}

dbl f20nl0_FIT_qqbar(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl f = -2.32235 * b * rho + 44.3522 * pow(b, 2) * rho - 24.6747 * pow(b, 3) * rho + 3.11129 * lrho * rho +
           2.92101 * b * pow(rho, 2) + 224.311 * pow(b, 2) * pow(rho, 2) + 21.5307 * pow(b, 3) * pow(rho, 2) +
           100.125 * lrho * pow(rho, 2) + 2.05531 * b * pow(rho, 3) + 945.506 * pow(b, 2) * pow(rho, 3) +
           36.1101 * pow(b, 3) * pow(rho, 3) - 176.632 * pow(b, 4) * pow(rho, 3) + 563.1 * lrho * pow(rho, 3) +
           7.68918 * pow(b, 4) * pow(rho, 4) + 568.023 * lrho * pow(rho, 4);
  return f;
}

//===========================================================================
//   Fit for the power suppressed terms (contains the qqbar "constant" term)
//   NOTE: the term ~nl^2 has no "threshold", and is known analytically
//===========================================================================

dbl PXS_NNLO_FIT_qqbar(cdbl rho, cdbl nl) {
  // double nl=5.0;
  // the functions fij are normalized in as^n/m^2
  cdbl f20 = pow(nl, 2) * f20nl2_FIT_qqbar(rho) + nl * f20nl1_FIT_qqbar(rho) + f20nl0_FIT_qqbar(rho);
  // double f21 = pow(nl,2)*f21nl2_FIT_qqbar(rho) + nl*f21nl1_FIT_qqbar(rho) + f21nl0_FIT_qqbar(rho);
  // double f22 = pow(nl,2)*f22nl2_FIT_qqbar(rho) + nl*f22nl1_FIT_qqbar(rho) + f22nl0_FIT_qqbar(rho);

  // double f = f20 + f21*(lg + lgfr) + f22*(pow(lg,2) + 2.*lg*lgfr + pow(lgfr,2));
  return f20;
}

dbl f2qqbar(cdbl rho, cdbl nl) {
  cdbl f = PXS_NNLO_THRESHOLD_qqbar(rho, nl) + PXS_NNLO_FIT_qqbar(rho, nl) + PXS_NNLO_identical_qqbar(rho);
  return f;
}

//=======================================================================
// The NNLO gq functions (i.e. the ones without scale-logs):
//=======================================================================

dbl f20nl0_FIT_gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lbe = log(b);
  cdbl lrho2 = pow(lrho, 2);
  cdbl lrho3 = pow(lrho, 3);
  cdbl lbe2 = pow(lbe, 2);
  cdbl rho2 = pow(rho, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);
  cdbl b7 = pow(b, 7);

  cdbl f = 28.0998 * b2 + 24.1753 * b3 - 12.3211 * b5 - 49.909 * b7 + 11.7853 * b3 * lbe + 28.6697 * b6 * lbe2 -
           3.82993 * lrho3 * rho + lrho2 * (-9.80339 * rho - 76.7407 * rho2) + lrho * (-1.68957 + 30.6335 * rho2);
  return f;
}

dbl f20nl1_FIT_gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lbe = log(b);
  cdbl lrho2 = pow(lrho, 2);
  cdbl lrho3 = pow(lrho, 3);
  cdbl lbe2 = pow(lbe, 2);
  cdbl rho2 = pow(rho, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b7 = pow(b, 7);

  cdbl f = 0.363838 * b2 - 1.44391 * b3 + 1.1146 * b7 - 0.309165 * b3 * lbe + 0.990057 * b4 * lbe2 +
           0.0401411 * lrho3 * rho + 0.362183 * lrho * rho2 + lrho2 * (0.194867 * rho + 1.57274 * rho2);
  return f;
}

dbl f2gq(cdbl rho, cdbl nl) {
  // double nl=5.0;
  // the functions fij are normalized as: as^n/m^2
  cdbl f20 = f20nl0_FIT_gq(rho) + nl * f20nl1_FIT_gq(rho);
  return f20;
}

//==========================================================================
// The scaling functions ~Log(muF^2/m^2) for gq-reaction.
//
// Implemented fits: as in ``Hathor" (Aliev et al. arXiv:1007.1327v1)
//==========================================================================

dbl f21nl0_FIT_gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lbe = log(b);
  cdbl lrho2 = pow(lrho, 2);
  cdbl lbe2 = pow(lbe, 2);
  cdbl b2 = pow(b, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);
  cdbl b7 = pow(b, 7);

  cdbl f =
      -0.00272708 * b2 - 0.194727 * b3 - 7.74780189742068692567 * b4 - 32.9814144136482376215 * b5 -
      2168.37139107242508102 * b6 + 2213.66358353965767595 * b7 + 0.62267243714534082136754351670760392811 * b3 * lbe -
      1.469367779446109471946 * b4 * lbe - 82.4683103579544036772 * b5 * lbe - 739.731611466006394139 * b6 * lbe -
      1201.795933886655196561 * b7 * lbe - 0.56735789898499726176870600831868545725 * b3 * lbe2 +
      215.848631557653764823 * b3 * lrho - 1146.873278298064968087 * b4 * lrho + 1744.26195133738780366 * b5 * lrho -
      1000.346862873385464604 * b6 * lrho + 187.933944366753454585 * b7 * lrho + 291.663849809244439603 * b3 * lrho2 -
      772.77269566517296654 * b4 * lrho2 + 653.225988602083032794 * b5 * lrho2 - 172.118699492821641978 * b6 * lrho2;
  return f;
}

dbl f21nl1_FIT_gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lbe = log(b);
  cdbl lrho2 = pow(lrho, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);
  cdbl b7 = pow(b, 7);

  cdbl f = 0.00554947127181299611557 * b3 + 0.0225430311833103251235 * b4 - 0.63483974233126651067 * b5 +
           9.99566887547610873574 * b6 - 9.40032480635199279045 * b7 -
           0.0112980360836839281564331066591599701444 * b3 * lbe + 0.00425431823252394833895 * b4 * lbe +
           0.252351741613262253035 * b5 * lbe + 2.40450791609598725401 * b6 * lbe + 4.53123563973969674337 * b7 * lbe -
           1.382674579044874428 * b3 * lrho + 6.3852062287204969368 * b4 * lrho - 8.30987940966125596394 * b5 * lrho +
           3.65509528520063341522 * b6 * lrho - 0.347592328515009762391 * b7 * lrho -
           1.263737185806346314577 * b3 * lrho2 + 3.10439060379970187549 * b4 * lrho2 -
           2.30142593570054254576 * b5 * lrho2 + 0.460779458031490856105 * b6 * lrho2;
  return f;
}

dbl f22nl0_FIT_gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lbe = log(b);
  cdbl lrho2 = pow(lrho, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);
  cdbl b7 = pow(b, 7);

  cdbl f = -0.0956665797118411055 * b3 + 0.0796230395754896639 * b4 - 4.65135359743332977 * b5 +
           54.2282042897508657 * b6 - 50.3664336690071789 * b7 + 0.094559649830832876961451001386447576209 * b3 * lbe +
           0.01454877297063897196 * b4 * lbe + 1.090238571394724594 * b5 * lbe + 12.39185329777043355 * b6 * lbe +
           25.3325922738649971 * b7 * lbe - 8.19094150408542609 * b3 * lrho + 34.1888813941803399 * b4 * lrho -
           44.3046767125567972 * b5 * lrho + 20.6517749275371715 * b6 * lrho - 2.50412783544030659 * b7 * lrho -
           6.09708063744189756 * b3 * lrho2 + 16.0087447100660202 * b4 * lrho2 - 12.77595728730197071 * b5 * lrho2 +
           2.8643356372167124 * b6 * lrho2;
  return f;
}

dbl f22nl1_FIT_gq(cdbl rho) {
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);
  cdbl lbe = log(b);
  cdbl lrho2 = pow(lrho, 2);
  cdbl b3 = pow(b, 3);
  cdbl b4 = pow(b, 4);
  cdbl b5 = pow(b, 5);
  cdbl b6 = pow(b, 6);
  cdbl b7 = pow(b, 7);

  cdbl f = 0.00282450625811508503 * b3 + 0.0001125529801773163916 * b4 - 0.0325005325377584271 * b5 +
           0.171745428117885934 * b6 - 0.1362549485229853436 * b7 + 0.0000160504261251611439 * b4 * lbe +
           0.00258676966103027447 * b5 * lbe + 0.0358936097435228275 * b6 * lbe + 0.0721640494225125838 * b7 * lbe -
           0.0404568351507943864 * b3 * lrho + 0.1116658960293489047 * b4 * lrho - 0.1145615544900474319 * b5 * lrho +
           0.0537148479428972514 * b6 * lrho - 0.01036242513667345341 * b7 * lrho - 0.0253066537137784997 * b3 * lrho2 +
           0.0551710547372787562 * b4 * lrho2 - 0.0403074509370066612 * b5 * lrho2 +
           0.01044304642981058897 * b6 * lrho2;
  return f;
}

dbl fbarF2gq(cdbl rho, cdbl nl) {
  cdbl f21 = f21nl0_FIT_gq(rho) + nl * f21nl1_FIT_gq(rho);
  return f21;
}

dbl fbarR2gq(cdbl rho, cdbl nl) {
  cdbl f = 3. * beta0(nl) * f1gq(rho, nl);
  return f;
}

dbl fbarRF2gq(cdbl rho, cdbl nl) {
  cdbl f = 3. * beta0(nl) * fbarF1gq(rho, nl);
  return f;
}

dbl fbarFF2gq(cdbl rho, cdbl nl) {
  cdbl f22 = f22nl0_FIT_gq(rho) + nl * f22nl1_FIT_gq(rho);
  return f22 - fbarRF2gq(rho, nl);
}

//=======================================================================
// The NNLO functions (no scale-logs) for qq-, qq'- and qqbar'-reactions
//=======================================================================

dbl f2qq(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -0.740558 - 22.8129 * pow(b, 2) - 0.191648 * pow(b, 3) - 6.58031 * pow(b, 4) - 0.537669 * pow(b, 5) +
              31.7872 * pow(b, 6) - 21.0783 * lrho * rho - 3.25313 * pow(lrho, 2) * rho - 10.8176 * lrho * pow(rho, 2) +
              15.8988 * pow(lrho, 2) * pow(rho, 2) + 8.64557 * lrho * pow(rho, 3);
  return -0.4768323995789214 * lrho - pow(b, 2) * exp(fexp);
}

dbl f2qqprime(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -0.740558 - 23.4518 * pow(b, 2) - 0.193073 * pow(b, 3) - 5.97215 * pow(b, 4) - 0.541402 * pow(b, 5) +
              31.8227 * pow(b, 6) - 21.3725 * lrho * rho - 3.29162 * pow(lrho, 2) * rho - 11.1642 * lrho * pow(rho, 2) +
              15.9932 * pow(lrho, 2) * pow(rho, 2) + 8.64746 * lrho * pow(rho, 3);
  return -0.4768323995789214 * lrho - pow(b, 2) * exp(fexp);
}

dbl f2qqbarprime(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -0.740572 - 31.2117 * pow(b, 2) - 0.31495 * pow(b, 3) + 15.8601 * pow(b, 4) - 1.64639 * pow(b, 5) +
              18.9767 * pow(b, 6) - 19.6977 * lrho * rho - 3.16565 * pow(lrho, 2) * rho - 16.1386 * lrho * pow(rho, 2) +
              12.3828 * pow(lrho, 2) * pow(rho, 2) + 4.17707 * lrho * pow(rho, 3);
  return -0.4768323995789214 * lrho - pow(b, 2) * exp(fexp);
}

//==========================================================================
// The scaling functions ~Log(muF^2/m^2) for qq-, qq'- and qqbar'-reactions
//==========================================================================

// same for all 3 reactions
dbl fbarFF2qqprime(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -2.20937 - 99.5512 * pow(b, 2) + 0.394393 * pow(b, 3) + 92.3169 * pow(b, 4) + 3.9297 * pow(b, 5) +
              4.81762 * pow(b, 8) - 43.9374 * lrho * rho - 6.31738 * pow(lrho, 2) * rho - 48.7898 * lrho * pow(rho, 2) +
              25.8588 * pow(lrho, 2) * pow(rho, 2) - 7.30711 * lrho * pow(rho, 3);
  cdbl f = -0.157185 * pow(b, 2) - 0.047419 * lrho + pow(b, 2) * pow(rho, 0.72854) * exp(fexp);
  return f;
}

// same for qq' and qqbar'
dbl fbarF2qqprime(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -0.07924228 + 211.8104 * pow(b, 2) - 0.4464898 * pow(b, 3) - 196.6873 * pow(b, 4) - 5.141548 * pow(b, 5) -
              9.448772 * pow(b, 8) + 90.67104 * lrho * rho + 11.02778 * pow(lrho, 2) * rho +
              104.2288 * lrho * pow(rho, 2) - 55.31094 * pow(lrho, 2) * pow(rho, 2) + 16.62298 * lrho * pow(rho, 3);
  cdbl f = 1.17852 * pow(b, 2) + 0.254683 * lrho - 1. * pow(b, 2) * pow(rho, 0.4) * exp(fexp);
  return f;
}

// for qq only
dbl fbarF2qq(cdbl rho, cdbl nl) {
  (void)nl;
  cdbl b = sqrt(1 - rho);
  cdbl lrho = log(rho);

  cdbl fexp = -0.07923334 + 209.7125 * pow(b, 2) - 0.3786914 * pow(b, 3) - 195.7483 * pow(b, 4) - 4.500831 * pow(b, 5) -
              9.000056 * pow(b, 8) + 89.56873 * lrho * rho + 10.92725 * pow(lrho, 2) * rho +
              103.6116 * lrho * pow(rho, 2) - 54.04722 * pow(lrho, 2) * pow(rho, 2) + 16.25028 * lrho * pow(rho, 3);
  cdbl f = 1.17852 * pow(b, 2) + 0.254683 * lrho - 1. * pow(b, 2) * pow(rho, 0.4) * exp(fexp);
  return f;
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
  /** @brief NNLO */
  coeff f2;
  /** @brief NNLO F SV */
  coeff fbarF2;
  /** @brief NNLO R SV */
  coeff fbarR2;
  /** @brief NNLO FF SV */
  coeff fbarFF2;
  /** @brief NNLO RF SV */
  coeff fbarRF2;
  /** @brief NNLO RR SV */
  coeff fbarRR2;
};

/** @brief gluon-gluon channel */
const CoeffMap gg = {f0gg, f1gg, fbarF1gg, fbarR1gg, f2gg, 0, 0, 0, 0, 0};

/** @brief quark-antiquark channel */
const CoeffMap qqbar = {f0qqbar, f1qqbar, fbarF1qqbar, fbarR1qqbar, f2qqbar, 0, 0, 0, 0, 0};

/** @brief gluon-quark channel */
const CoeffMap gq = {0, f1gq, fbarF1gq, 0, f2gq, fbarF2gq, fbarR2gq, fbarFF2gq, fbarRF2gq, 0};

/** @brief quark-quark channel */
const CoeffMap qq = {0, 0, 0, 0, f2qq, fbarF2qq, 0, fbarFF2qqprime, 0, 0};

/** @brief quark-antiquark' channel */
const CoeffMap qqbarprime = {0, 0, 0, 0, f2qqbarprime, fbarF2qqprime, 0, fbarFF2qqprime, 0, 0};

/** @brief quark-quark' channel */
const CoeffMap qqprime = {0, 0, 0, 0, f2qqprime, fbarF2qqprime, 0, fbarFF2qqprime, 0, 0};

///@}

}  // namespace FixedOrder

}  // namespace MaunaKea

#endif  // FO_HPP_
