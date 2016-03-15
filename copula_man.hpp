#ifndef COPULA_MAN_HPP
#define COPULA_MAN_HPP

/*---------------------------------------------*/
/* Bivariate copula manipulation               */
/* Lei Hua, Mar, 2016                          */
/*---------------------------------------------*/

typedef double (* Cfun)(double u, double v, std::vector<double> par);
double dMixKT(double u, double v, std::vector<double> par, Cfun pC, Cfun C1, Cfun C2, Cfun dC);
double C1MixKT(double u, double v, std::vector<double> par, Cfun pC, Cfun C1);
double C2MixKT(double u, double v, std::vector<double> par, Cfun pC, Cfun C1);
double invC2B6(double u, double v, std::vector<double> par);

#endif /* COPULA_MAN_HPP */

