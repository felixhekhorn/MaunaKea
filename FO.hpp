#ifndef FO_HPP_
#define FO_HPP_

#include "./config.h"

namespace MaunaKea {

double f0gg(cdbl rho) {
  cdbl beta = sqrt(1 - rho);
  cdbl f = M_PI * beta * (rho / 192.) *
           ((pow(rho, 2) + 16 * rho + 16) * log((1 + beta) / (1 - beta)) / beta - 28 - 31 * rho);
  return f;
}

double f0qqbar(cdbl rho) {
  cdbl beta = sqrt(1 - rho);
  return M_PI * beta * rho / 27. * (2 + rho);
}

}  // namespace MaunaKea

#endif  // FO_HPP_
