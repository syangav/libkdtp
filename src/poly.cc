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
#include <cmath>
#include <vector>

#include "kdtp.h"

static std::vector<double> cardano(double p, double q);
static std::vector<double> ferrari(double p, double q, double r);


/**
 * Computing real roots of a second order polynomial.
 *
 * Input
 *    double a : coefficient of the second order
 *    double b : coefficient of the first order
 *    double c : constant term
 *
 * Output
 *    std::vector<double> : real solutions a*X^2+b*X+c=0
 */
std::vector<double>
kdtp::poly_root_2(double a, double b, double c)
{
  std::vector<double> sol;

  if (fabs(a) < kdtp::EPSILON) {
    if (fabs(b) < kdtp::EPSILON)
      sol.push_back(-c);
    else
      sol.push_back(-c/b);
  } else {
    double d = b*b-4*a*c;

    if (fabs(d) < kdtp::EPSILON)
      sol.push_back(-b/(2*a));
    else if (d > 0) {
      sol.push_back((-b-sqrt(d))/(2*a));
      sol.push_back((-b+sqrt(d))/(2*a));
    }
  }
  return sol;
}


/**
 * Computing real roots of a third order polynomial using Cardano's method.
 *
 * Input
 *    double a : coefficient of the third order
 *    double b : coefficient of the second order
 *    double c : coefficient of the first order
 *    double d : constant term
 *
 * Output
 *    std::vector<double> : real solutions of a*X^3+b*X^2+c*X+d=0
 */
std::vector<double>
kdtp::poly_root_3(double a, double b, double c, double d)
{
  if (fabs(a) < kdtp::EPSILON) return poly_root_2(b, c, d);

  double aa = b/a;
  double bb = c/a;
  double cc = d/a;

  double p = bb-aa*aa/3.;
  double q = 2*aa*aa*aa/27.-aa*bb/3.+cc;
  std::vector<double> sol = cardano(p, q);

  for(unsigned int k = 0; k<sol.size(); k++)
    sol[k] = sol[k]-aa/3;

  return sol;
}


/**
 * Computing real roots of a fourth order polynomial using Ferrari's method.
 *
 * Input
 *    double a : coefficient of the fourth order
 *    double b : coefficient of the third order
 *    double c : coefficient of the second order
 *    double d : coefficient of the first order
 *    double e : constant term
 *
 * Output
 *    std::vector<double> : real solutions of a*X^4+b*X^3+c*X^2+d*X+e=0
 */
std::vector<double>
kdtp::poly_root_4(double a, double b, double c, double d, double e)
{
  if (fabs(a) < kdtp::EPSILON) return poly_root_3(b, c, d, e);

  double aa = b/a;
  double bb = c/a;
  double cc = d/a;
  double dd = e/a;
  double p = bb-3*aa*aa/8;
  double q = cc-aa*bb/2+aa*aa*aa/8;
  double r = dd-aa*cc/4+bb*aa*aa/16-3*aa*aa*aa*aa/256;

  std::vector<double> sol = ferrari(p, q, r);
  for(unsigned int k = 0; k<sol.size(); k++)
    sol[k] = sol[k] - aa/4.;

  return sol;
}


/**
 * Implementation of Cardano's method.
 *
 * Input :
 *    double p : coefficient of the first order
 *    double q : constant term
 *
 * Output :
 *    std::vector<double> : real solutions of X^3+p*X+q=0
 */
static std::vector<double>
cardano(double p, double q)
{
  std::vector<double> sol;

  double delta = -(4*p*p*p+27*q*q);
  if (fabs(delta) < kdtp::EPSILON) {
    if (fabs(p)<kdtp::EPSILON || fabs(q) < kdtp::EPSILON)
      sol.push_back(0);
    else {
      sol.push_back(3*q/p);
      sol.push_back(-3*q/(2*p));
    }
  } else if (delta < 0) {
    double v = sqrt(-delta/27)/2;
    double u = -q/2+v;

    u = kdtp::sign(u)*pow(fabs(u),1/3.);
    v = -q/2-v;
    v = kdtp::sign(v)*pow(fabs(v),1/3.);
    sol.push_back(u+v);
  } else {
    double u = 1/3. * acos((-q/2)*sqrt(27/(-p*p*p)));
    double v=2*sqrt(-p/3);
    for(int k = 0; k<3; k++)
      sol.push_back(v*cos(u+2*((double)k)*M_PI/3));
  }

  return sol;
}

/**
 * Implementation of Ferrari's method.
 *
 * Input
 *    double p : coefficient of the second order
 *    double q : coefficient of the first order
 *    double r : constant term
 *
 * Output
 *    std::vector<double> : real solutions of X^4+p*X^2+q*X+r=0
 */
static std::vector<double>
ferrari(double p, double q, double r)
{
  std::vector<double> sol;
  if (fabs(q) < kdtp::EPSILON) {
    std::vector<double> soltmp = kdtp::poly_root_2(1, p, r);
    for(unsigned k = 0; k<soltmp.size(); k++) {
      if (soltmp.at(k) >= 0.) {
        sol.push_back(-sqrt(soltmp.at(k)));
        sol.push_back( sqrt(soltmp.at(k)));
      }
    }
  } else {
    std::vector<double> soltmp = kdtp::poly_root_3(1, -p/2, -r, p*r/2-q*q/8);
    double y0 = soltmp.at(0);
    double delta = 2*y0-p;
    if (delta < 0.)
      sol = kdtp::poly_root_2(1, 0, y0);
    else {
      double a0 = sqrt(delta);
      double b0;
      if (a0 < kdtp::EPSILON) {
        if (y0*y0-r>0)
          b0 = sqrt(y0*y0-r);
        else
          b0 = 0.;
      } else
        b0 = -q/(2*a0);

      sol = kdtp::poly_root_2(1, a0, y0+b0);
      soltmp = kdtp::poly_root_2(1, -a0, y0-b0);
      for(unsigned int k = 0; k<soltmp.size(); k++)
        sol.push_back(soltmp.at(k));
    }
  }

  return sol;
}
