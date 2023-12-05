/**
 * @file config.h
 * @brief includes all common header files and defines types and shortcuts
 */

#ifndef Config_H_
#define Config_H_

#include <stdlib.h>
#include <cmath>
#include <string>
#include <stdexcept>
#include <iostream>
#include <stdio.h>

using namespace std;

/** @name shorthands */
///@{
/** @brief unsigned int shorthand */
typedef unsigned int uint;
/** @brief const unsigned int shorthand */
typedef const uint cuint;
/** @brief const int shorthand */
typedef const int cint;

/** @brief string shorthand */
typedef string str;

/** @brief defines floating point precision */
typedef double dbl;
/** @brief defines floating point precision */
typedef const dbl cdbl;

/** @brief NaN (often (mis-)used by me as synonym to null-pointer/void) */
cdbl dblNaN = nan("");
///@}

#endif // Config_H_
