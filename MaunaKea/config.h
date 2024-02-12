/**
 * @brief includes all common header files and defines types and shortcuts
 */

#ifndef MAUNAKEA_CONFIG_H_
#define MAUNAKEA_CONFIG_H_

#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

/** @name shorthands */
///@{
/** @brief unsigned int shorthand */
typedef unsigned int uint;
/** @brief const unsigned int shorthand */
typedef const uint cuint;
/** @brief const int shorthand */
typedef const int cint;

/** @brief string shorthand */
typedef std::string str;

/** @brief defines floating point precision */
typedef double dbl;
/** @brief defines floating point precision */
typedef const dbl cdbl;

/** @brief NaN (often (mis-)used by me as synonym to null-pointer/void) */
cdbl dblNaN = std::nan("");
///@}

#endif  // MAUNAKEA_CONFIG_H_
