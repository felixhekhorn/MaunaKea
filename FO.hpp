#include "config.h"

double f0gg(cdbl rho)
{
  cdbl beta = sqrt(1-rho);
  cdbl f = M_PI*beta*(rho/192.)
    *((pow(rho,2) + 16*rho+16)
      *log((1+beta)/(1-beta))/beta-28-31*rho);	
  return f;
};
