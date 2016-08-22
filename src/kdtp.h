/*
 * Copyright (c) 2016 LAAS/CNRS
 * All rights reserved.
 *
 * Redistribution and use  in source  and binary  forms,  with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *
 *   1. Redistributions of  source  code must retain the  above copyright
 *      notice and this list of conditions.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice and  this list of  conditions in the  documentation and/or
 *      other materials provided with the distribution.
 *
 *                                     Alexandre Boeuf on Fri Aug 19 2016
 */

#ifndef H_KDTP
#define H_KDTP

#include <cmath>

#include "libkdtp.h"

namespace kdtp {

  static const double EPSILON = 1e-10;
  static const double EPSI4 = 1e-4;
  static const unsigned int NB_IT_MAX = 1000;

  static inline int sign(double d)
  {
    return (d > EPSILON) - (d < -EPSILON);
  }

  static inline int sqr(double d)
  {
    return d * d;
  }

  std::vector<double> poly_root_2(double a, double b, double c);
  std::vector<double> poly_root_3(double a, double b, double c, double d);
  std::vector<double> poly_root_4(double a, double b, double c, double d,
                                  double e);
}


#endif /* H_KDTP */
