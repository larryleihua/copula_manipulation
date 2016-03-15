/*---------------------------------------------*/
/* Copula manipulation                         */
/* Lei Hua, March 2016                         */
/*---------------------------------------------*/

#include <iostream>
#include <functional>
#include <vector>
#include <numeric>
#include <algorithm>
#include <Rcpp.h>
#include "copula_man.hpp"

#ifdef STUDENT_T
#include <boost/math/distributions/students_t.hpp>
#endif

#define UMAX  1-1e-10
#define UMIN  1e-10
#define XINFMAX DBL_MAX
#define TOL 1e-6
#define MAXIT 1000


// [[Rcpp::depends(RcppGSL)]]

// [[Rcpp::depends(RcppArmadillo)]]

using std::pow; using std::log;
using std::exp; using std::vector;


// error handling for boost::math

#ifdef STUDENT_T
using boost::math::policies::policy;
using boost::math::policies::overflow_error;
using boost::math::policies::underflow_error;
using boost::math::policies::domain_error;
using boost::math::policies::pole_error;
using boost::math::policies::denorm_error;
using boost::math::policies::evaluation_error;
using boost::math::policies::ignore_error;

// Define a custom policy to ignore just overflow:
typedef policy< overflow_error<ignore_error> > ignore_overflow;

// ignore all errors: domain, pole, overflow, underflow, denorm & evaluation:
typedef policy< domain_error<ignore_error>,
                pole_error<ignore_error>,
				overflow_error<ignore_error>,
				underflow_error<ignore_error>,
				denorm_error<ignore_error>,
				evaluation_error<ignore_error>
				> ignore_allerror;

typedef boost::math::students_t_distribution<double, ignore_allerror> my_students_t;
#endif
				
////////////////////////////////////////
//     copula parameters ranges       //
////////////////////////////////////////

// Student t: mu>0, -1 < r < 1
// B4:  de>0
// B5: 	de>1
// B6:  de>1
// BB1: de>=1, th>0
// BB7:  de>0, th>=1

// initializing the copula parameters based on the range of them
vector<double> cparIni(vector<double> cpar, const int family)
{
    
    vector<double> out;
    switch( family )
    {
        case 4:  // B4
            out.push_back(exp(cpar[0]));
            break;
        case 5:  // B5
            out.push_back(exp(cpar[0])+1);
            break;
        case 6:  // B6
            out.push_back(exp(cpar[0])+1);
            break;
        default:
            break;
    }
    return out;
}


//////////////////////////////////////////
// Standard Bivariate copulas           //
//--------------------------------------//
// pC: cdf of C                         //
// C1: derivative to u1 of C            //
// C2: derivative to u2 of C            //
// dC: density of C                     // 
//////////////////////////////////////////

#ifdef STUDENT_T
// Bivariate Student t copula
double dSt(double u, double v, vector<double> par)
{
  double nu = par[0]; // here nu can be non-integral
  double r = par[1];
  double out = 0;
  
  double rho2,h1,h2,h3,h4,h5,h6;
  double x,y,x2,y2;
  
  rho2 = pow(r,2);
  h1 = 1-rho2;
  h2 = nu/2;
  h3 = h2 + 0.5;
  h4 = h2 + 1;
  h5 = 1/ nu;
  h6 = h5/h1;
  
  my_students_t dist(nu);
  
  x = boost::math::quantile(dist, u);
  y = boost::math::quantile(dist, v);
  x2 = pow(x,2);
  y2 = pow(y,2);
  
  out = tgamma(h4)*tgamma(h2)/sqrt(h1)/pow(tgamma(h3),2)*pow(1+h5*x2,h3)*pow(1+h5*y2,h3)/pow(1+h6*(x2+y2-2*r*x*y),h4);
  return out;
}
#endif

// based on CopulaModel R package
// xq: nodes;  wq: weights
#define EPS 3.0e-11
void gauleg(int nq, vector<double>& xq, vector<double>& wq)
{
    size_t m,j,i,n;
    double z1,z,xm,xl,pp,p3,p2,p1;
    n = nq; 
    m = (n+1)/2;
  
    // boundary points 
    double x1 = 0;
    double x2 = 1;
    
    xm = 0.5*(x2 + x1);
    xl = 0.5*(x2 - x1);
    for(i=1;i<=m;++i) // yes, i starts from 1
    {
        z = cos(3.14159265358979323846*(i-0.25)/(n+0.5));
        do
        {
            p1=1.0;
            p2=0.0;
            for(j=1;j<=n;j++)
            {
                p3=p2;
                p2=p1;
                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
            }
            pp=n*(z*p1-p2)/(z*z-1.0);
            z1=z;
            z=z1-p1/pp;
        }while (fabs(z-z1) > EPS);
        
        xq[i-1] = xm-xl*z;
        xq[n-i]=xm+xl*z;
        wq[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
        wq[n-i]=wq[i-1];
        
  }
}
#undef EPS

#ifdef STUDENT_T
// for cdf, open sources only have integral nu
// here we implement using numerical integration
double pSt(double u, double v, vector<double> par)
{
    double out = 0;
    int i, j;
    int nq = 27;

    vector<double> xl(nq), wl(nq);
    gauleg(nq, xl, wl);

    Rcpp::NumericMatrix intg(nq, nq);
    
    for(i=0; i < nq; ++i)
    {
        for(j=0; j < nq; ++j)
        {
            intg(i,j) = dSt(xl[i]*u, xl[j]*v, par)*u*v;
            out += wl[i]*wl[j]*intg(i,j);
        }
    }
    return out;
}

// integrand used to obtain C1St
// to-do: needs to be corrected
double _C1St_intg(double u, double v, vector<double> par)
{
    double nu = par[0];
    double r = par[1];
    double tem1, tem2, x1, x2, xR1x, out;
    double r2 = pow(r,2);
    
    boost::math::students_t dist(nu);
    x1 = boost::math::quantile(dist, u);
    x2 = boost::math::quantile(dist, v);
    
    xR1x = x1*(-x1/(-1+r2)+x2*r/(-1+r2))+x2*(x1*r/(-1+r2)-x2/(-1+r2));
    
    tem1 = tgamma((nu+2)/2) * tgamma(nu/2) / pow(tgamma((nu+1)/2),2) / pow(1-r2, 0.5);
    tem2 = pow((1 + pow(x1,2)/nu ), 1+nu) / pow(1+ xR1x / nu,(nu+2)/2);
    out = tem1 * tem2 * boost::math::pdf(dist, x1);  
    return out;
}

double C1St_v2(double u, double v, vector<double> par)
{
    double out = 0;
    
    int i;
    int nq = 27;

    vector<double> xl(nq), wl(nq);
    gauleg(nq, xl, wl);

    for(i=0; i < nq; ++i)
    {
        out +=  wl[i] * _C1St_intg(u, xl[i]*v, par)*v;
    }
    return out;
}


double C1St(double u, double v, vector<double> par)
{
    double out = 0;
    
    int i;
    int nq = 27;

    vector<double> xl(nq), wl(nq);
    gauleg(nq, xl, wl);

    for(i=0; i < nq; ++i)
    {
        out +=  wl[i] * dSt(u, xl[i]*v, par)*v;
    }
    return out;
}

double C2St(double u, double v, vector<double> par)
{
    double out = C1St(v, u, par);
    return out;
}
#endif


// B4 copula
double pB4(double u, double v, vector<double> par)
{
    double t1, t2, out;
    double de = par[0];
    t1 = pow(u,-1.0*de);
    t2 = pow(v,-1.0*de);
    out = pow(t1+t2-1.0,-1.0/de);
    return out;
}

double C1B4(double u, double v, vector<double> par)
{
    double t1, t2, t3, t5, out;
    double de = par[0];
    t1 = pow(u,-1.0*de);
    t2 = pow(v,-1.0*de);
    t3 = t1+t2-1.0;
    t5 = pow(t3,-1.0/de);
    out = t5*t1/u/t3;
    return out;
}

double C2B4(double u, double v, vector<double> par)
{
    double t1, t2, t3, t5, out;
    double de = par[0];
    t1 = pow(u,-1.0*de);
    t2 = pow(v,-1.0*de);
    t3 = t1+t2-1.0;
    t5 = pow(t3,-1.0/de);
    out = t5*t2/v/t3;
    return out;
}

double dB4(double u, double v, vector<double> par)
{
    double t1, t2, t3, t5, t7, t9, t10, t12, out;
    double de = par[0];
    t1 = pow(u,-1.0*de);
    t2 = pow(v,-1.0*de);
    t3 = t1+t2-1.0;
    t5 = pow(t3,-1.0/de);
    t7 = 1/v;
    t9 = t3*t3;
    t10 = 1/t9;
    t12 = 1/u;
    out = t5*t2*t7*t10*t1*t12+t5*t1*t12*t10*t2*de*t7;
    return out;
}

// B5 copula
double pB5(double u, double v, vector<double> par)
{
    double t2, t4, t8, out;
    double de = par[0];
    t2 = pow(1.0-u,1.0*de);
    t4 = pow(1.0-v,1.0*de);
    t8 = pow(t2+t4-t2*t4,1.0/de);
    out = 1.0-t8;
    return out;
}

double C1B5(double u, double v, vector<double> par)
{
    double t1, t2, t4, t6, t7, t8, t10, t11, out;
    double de = par[0];
    t1 = 1.0-u;
    t2 = pow(t1,1.0*de);
    t4 = pow(1.0-v,1.0*de);
    t6 = t2+t4-t2*t4;
    t7 = 1/de;
    t8 = pow(t6,1.0*t7);
    t10 = t2*de;
    t11 = 1/t1;
    out = -t8*t7*(-t10*t11+t10*t11*t4)/t6;
    return out;
}

double C2B5(double u, double v, vector<double> par)
{
    double t2, t3, t4, t5, t6, t7, t8, t11, out;
    double de = par[0];
    t2 = pow(1.0-u,1.0*de);
    t3 = 1.0-v;
    t4 = pow(t3,1.0*de);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t7 = 1/de;
    t8 = pow(t6,1.0*t7);
    t11 = 1/t3;
    out = -t8*t7*(-t4*de*t11+t5*de*t11)/t6;
    return out;
}

double dB5(double u, double v, vector<double> par)  // 1 <= de < oo
{
    double t1, t2, t3, t4, t5, t6, t7, t8, t9, t13, t18, t21, t22, t24, t27, out;
    double de = par[0];
    t1 = 1.0-u;
    t2 = pow(t1,1.0*de);
    t3 = 1.0-v;
    t4 = pow(t3,1.0*de);
    t5 = t2*t4;
    t6 = t2+t4-t5;
    t7 = 1/de;
    t8 = pow(t6,1.0*t7);
    t9 = de*de;
    t13 = 1/t3;
    t18 = t6*t6;
    t21 = t2*de;
    t22 = 1/t1;
    t24 = t22*t4;
    t27 = (-t4*de*t13+t5*de*t13)/t18*(-t21*t22+t21*t24);
    out = -t8/t9*t27+t8*de*t2*t24*t13/t6+t8*t7*t27;
    return out;
}

// B6 copula
double pB6(double u, double v, vector<double> par)
{
    double t1, t2, t3, t4, t7, out;
    double de = par[0];
    t1 = log(u);
    t2 = pow(-t1,1.0*de);
    t3 = log(v);
    t4 = pow(-t3,1.0*de);
    t7 = pow(t2+t4,1.0/de);
    out = exp(-t7);
    return out;
}

double C1B6(double u, double v, vector<double> par)
{
    double t1, t2, t3, t4, t5, t7, t14, out;
    double de = par[0];
    t1 = log(u);
    t2 = pow(-t1,1.0*de);
    t3 = log(v);
    t4 = pow(-t3,1.0*de);
    t5 = t2+t4;
    t7 = pow(t5,1.0/de);
    t14 = exp(-t7);
    out = -t7*t2/u/t1/t5*t14;
    return out;
}

double C2B6(double u, double v, vector<double> par)
{
    double t1, t2, t3, t4, t5, t7, t14, out;
    double de = par[0];
    t1 = log(u);
    t2 = pow(-t1,1.0*de);
    t3 = log(v);
    t4 = pow(-t3,1.0*de);
    t5 = t2+t4;
    t7 = pow(t5,1.0/de);
    t14 = exp(-t7);
    out = -t7*t4/v/t3/t5*t14;
    return out;
}

double dB6(double u, double v, vector<double> par)  // 1 <= de < oo
{
    double t1, t2, t3, t4, t5, t7, t9, t10, t11, t13, t14, t18, t19, t31, out;
    double de = par[0];
    t1 = log(u);
    t2 = pow(-t1,1.0*de);
    t3 = log(v);
    t4 = pow(-t3,1.0*de);
    t5 = t2+t4;
    t7 = pow(t5,1.0/de);
    t9 = 1/v;
    t10 = 1/t3;
    t11 = t9*t10;
    t13 = t5*t5;
    t14 = 1/t13;
    t18 = 1/u/t1;
    t19 = exp(-t7);
    t31 = t7*t7;
    out = -t7*t4*t11*t14*t2*t18*t19+t7*t2*t18*t14*t19*t4*de*t9*t10+t31*t2*t18*t14*t4*t11*t19;
    return out;
}

// BB1 copula, see CopulaModel package
// [[Rcpp::export]]
double dBB1(double u, double v, vector<double> par)
{
    double de, th, de1, th1, ut, vt, x, y, sm, smd, tem, out;
    de = par[0];
    th = par[1];
    de1 = 1/de;
    th1 = 1/th;
    ut = pow(u,(-th))-1; 
    vt = pow(v,(-th))-1;
    x = pow(ut,de); 
    y = pow(vt,de);
    sm = x+y; 
    smd = pow(sm, de1);
    tem = pow((1+smd),(-th1-2)) * (th*(de-1)+(th*de+1)*smd);
    out = tem*smd*x*y*(ut+1)*(vt+1)/sm/sm/ut/vt/u/v;
    return out;
}

// BB7 copula, see CopulaModel package
// [[Rcpp::export]]
double dBB7(double u, double v, vector<double> par)
{
    double de, th, de1, th1, ut, vt, x, y, sm, smd,tem, out;
    de = par[0];
    th = par[1];
    de1 = 1/de;
    th1 = 1/th;
    ut = 1. - pow((1.-u), th);
    vt = 1. - pow((1.-v), th); 
    x = pow(ut, (-de)) - 1;
    y = pow(vt, (-de)) - 1;
    sm = x+y+1;
    smd = pow(sm, (-de1));
    tem = pow((1.-smd), (th1-2.));
    tem = tem*(th*(de+1)-(th*de+1)*smd);
    out = tem*smd*(x+1)*(y+1)*(1-ut)*(1-vt)/sm/sm/ut/vt/(1-u)/(1-v);
    return out;
}

double pBB1(double u, double v, vector<double> par)
{
    double tem1, tem2, tem3, out;
    double de = par[0];
    double th = par[1];
    
    tem1 = pow(pow(u, -th) - 1.0, de);
    tem2 = pow(pow(v, -th) - 1.0, de);
    tem3 = 1+pow((tem1 + tem2), 1/de);
    out = pow(tem3, -1 / th);
    return out;
}

double pBB7(double u, double v, vector<double> par)
{
    double t2, t4, t6, t8, t11, t14, out;
    double de = par[0];
    double th = par[1];
    
    t2 = pow(1.0-u,1.0*th);
    t4 = pow(1.0-t2,-1.0*de);
    t6 = pow(1.0-v,1.0*th);
    t8 = pow(1.0-t6,-1.0*de);
    t11 = pow(t4+t8-1.0,-1.0/de);
    t14 = pow(1.0-t11,1.0/th);
    out = 1.0-t14;
    return out;
}


// modified from R package CopulaModel
// C_{1|2}(u|v) for BB1 : based on derivative of transform variables
// par = copula parameter with (th,de) : th>0, de>1
// [[Rcpp::export]]
double C2BB1(double u, double v, vector<double> par)
{
    double de, th, de1, th1, vt, ut, x, y, sm, smd, tem, out;
    de = par[0];
    th = par[1];
    de1 = 1/de;
    th1 = 1/th;
    vt = pow(v, (-th)) - 1;
    ut = pow(u, (-th)) - 1;
    x = pow(vt, de);
    y = pow(ut, de);
    sm = x+y;
    smd = pow(sm, de1);
    tem = pow((1+smd), (-th1-1));
    out = tem*smd*x*(vt+1)/sm/vt/v;
    return out;
}

double C1BB1(double u, double v, vector<double> par)
{
    return C2BB1(v, u, par);
}


double C2BB7(double u, double v, vector<double> par)
{ 
    double de, th, de1, th1, vt, ut, x, y, sm, smd, tem, out;
    de = par[0];
    th = par[1];
    de1 = 1./de;
    th1 = 1./th;
    vt = 1.-pow((1.-v),th);
    ut = 1.-pow((1.-u),th); 
    x = pow(vt, -de)-1;
    y = pow(ut, -de)-1;
    sm = x+y+1;
    smd = pow(sm, -de1);
    tem = pow((1.-smd), (th1-1.));
    out = tem*smd*(x+1)*(1.-vt)/sm/vt/(1.-v);
    return out;
}

/*
double invC2BB1(double u, double v, vector<double> par)
{

    return out;
}
*/

//////////////////////////////////////////
// reflection-transform
// These functions are used to reflect a 
// bivariate copula through the line y = 1 - x
// Chat(u,v) = C(1-v, 1-u) + u + v - 1
// -------------------------------------------
// pRT:     cdf of reflection-transformed copula
// C1RT:    derivative to u1 of reflection-transformed copula
// C2RT:    derivative to u2 of reflection-transformed copula
// dRT:     density of reflection-transformed copula
//////////////////////////////////////////
double pRT(double u, double v, vector<double> par, Cfun pC)
{
    double out = pC(1-v, 1-u, par) + u + v - 1;
    return out;
}

// note that here C2 is used !!
double C1RT(double u, double v, vector<double> par, Cfun C2)
{
    double out = 1 - C2(1-v, 1-u, par);
    return out;
}

double C2RT(double u, double v, vector<double> par, Cfun C1)
{
    double out = 1 - C1(1-v, 1-u, par);
    return out;
}

double dRT(double u, double v, vector<double> par, Cfun dC)
{
    double out = dC(1-v, 1-u, par);
    return out;
}

// RT-B4:
double pRB4(double u, double v, vector<double> par)
{
    double out = pRT(u,v, par, pB4);
    return out;
}

double C1RB4(double u, double v, vector<double> par)
{
    double out = C1RT(u,v, par, C2B4);
    return out;
}

double C2RB4(double u, double v, vector<double> par)
{
    double out = C2RT(u,v, par, C1B4);
    return out;
}

// [[Rcpp::export]]
double dRB4(double u, double v, vector<double> par)
{
    double out = dRT(u,v, par, dB4);
    return out;
}

// RT-B5:
double pRB5(double u, double v, vector<double> par)
{
    double out = pRT(u,v, par, pB5);
    return out;
}

double C1RB5(double u, double v, vector<double> par)
{
    double out = C1RT(u,v, par, C2B5);
    return out;
}

double C2RB5(double u, double v, vector<double> par)
{
    double out = C2RT(u,v, par, C1B5);
    return out;
}

// [[Rcpp::export]]
double dRB5(double u, double v, vector<double> par)
{
    double out = dRT(u,v, par, dB5);
    return out;
}

// RT-B6:
double pRB6(double u, double v, vector<double> par)
{
    double out = pRT(u,v, par, pB6);
    return out;
}

double C1RB6(double u, double v, vector<double> par)
{
    double out = C1RT(u,v, par, C2B6);
    return out;
}

double C2RB6(double u, double v, vector<double> par)
{
    double out = C2RT(u,v, par, C1B6);
    return out;
}

// [[Rcpp::export]]
double dRB6(double u, double v, vector<double> par)
{
    double out = dRT(u,v, par, dB6);
    return out;
}

// RT-BB1:
double pRBB1(double u, double v, vector<double> par)
{
    double out = pRT(u,v, par, pBB1);
    return out;
}

double C1RBB1(double u, double v, vector<double> par)
{
    double out = C1RT(u,v, par, C2BB1);
    return out;
}

double C2RBB1(double u, double v, vector<double> par)
{
    double out = C2RT(u,v, par, C1BB1);
    return out;
}

double dRBB1(double u, double v, vector<double> par)
{
    double out = dRT(u,v, par, dBB1);
    return out;
}

//////////////////////////////////////////
// Khoudraji-transform
// --------------------------------------
// pKT: cdf of K-transformed bivariate copula
// dKT: density of K-transformed bivariate copula
// C1KT: derivative to u1 of pKT
// C2KT: derivative to u2 of pKT
// dMixKT: density of mixed K-transformed bivariate copula
// C1MixKT: derivative to u1 of mixed K-transformed bivariate copula
// C2MixKT: derivative to u2 of mixed K-transformed bivariate copula
//////////////////////////////////////////
double pKT(double u, double v, vector<double> par, Cfun pC)
{
    // par:     parameters for the original copula, and a1 a2 for KT
    // pC:      cdf of the original copula
 
    vector<double>::size_type nparC = par.size() - 2;     // nparC:   number of parameters of the original copula
    vector<double> parC(par.begin(), par.begin()+nparC );
    double a1 = par[nparC];
    double a2 = par[nparC+1];
    double out = pow(u, 1-a1) * pow(v, 1-a2) * pC( pow(u, a1), pow(v, a2), parC);
    return out;
}


double C1KT(double u, double v, vector<double> par, Cfun pC, Cfun C1)
{
    // par:     parameters for the original copula, and a1 a2 for KT
    // pC:      cdf of the original copula
    // C1:      derivative to u1 of the original copula
 
    vector<double>::size_type nparC = par.size() - 2;     // nparC:   number of parameters of the original copula
    vector<double> parC(par.begin(), par.begin()+nparC );
    double a1 = par[nparC];
    double a2 = par[nparC+1];
    double ua1 = pow(u, a1);
    double va2 = pow(v, a2);
    double tem1 = (1-a1) * pow(u, -a1) * pow(v, 1-a2) * pC(ua1, va2, parC);
    double tem2 = a1 * pow(v, 1-a2) * C1(ua1, va2, parC);
    double out = tem1 + tem2;
    return out;
}


double C2KT(double u, double v, vector<double> par, Cfun pC, Cfun C2)
{
    // par:     parameters for the original copula, and a1 a2 for KT
    // pC:      cdf of the original copula
    // C2:      derivative to u2 of the original copula
 
    vector<double>::size_type nparC = par.size() - 2;     // nparC:   number of parameters of the original copula
    vector<double> parC(par.begin(), par.begin()+nparC );
    double a1 = par[nparC];
    double a2 = par[nparC+1];
    double ua1 = pow(u, a1);
    double va2 = pow(v, a2);
    double tem1 = (1-a2) * pow(u, 1-a1) * pow(v, -a2) * pC(ua1, va2, parC);
    double tem2 = a2 * pow(u, 1-a1) * C2(ua1, va2, parC);
    double out = tem1 + tem2;
    return out;
}

double invC1MixKT(double u, double v, vector<double> par, Cfun pC, Cfun C1, Cfun C2, Cfun dC)
{
    const double tol = TOL;
     
    double CDF = v + 2*tol;
    double CDF1 = 100;
	double diff = 0.1;
    double lower, upper;
    double t;
    
    int maxiter = MAXIT;
    int kount = 0;
 
    /*--------------------------------
    # Now use modified Newton-Raphson
    #--------------------------------*/

    lower = 0;
    upper = 1;
    t = 0.5;
    
    while ( kount < maxiter && std::abs(v - CDF) > tol )
    {
        kount += 1;
		if( isnan(CDF) || isnan(CDF1) )
		{
			diff /=-2;
		}
		else
		{
			diff = (CDF-v) / CDF1;
		}	
        t = t - diff; 
        
        if (t < lower || t > upper)
        {
            t = 0.5 * (lower + upper);
        }

        CDF1 = dMixKT(u, t, par, pC, C1, C2, dC);
        CDF = C1MixKT(u, t, par, pC, C1);

        if ( CDF > v )
        {
            upper = t;
        }
        else
        {
            lower = t;
        }
    }
    return t;
}


double invC2MixKT(double u, double v, vector<double> par, Cfun pC, Cfun C1, Cfun C2, Cfun dC)
{
    const double tol = TOL;
     
    double CDF = v + 2*tol;
    double CDF1 = 100;
    double diff = 0.1;
    double lower, upper;
    double t;
    
    int maxiter = MAXIT;
    int kount = 0;
 
    /*--------------------------------
    # Now use modified Newton-Raphson
    #--------------------------------*/

    lower = 0;
    upper = 1;
    t = 0.5;
    
    while ( kount < maxiter && std::abs(u - CDF) > tol )
    {
        kount += 1;
		if( isnan(CDF) || isnan(CDF1) )
		{
			diff /=-2;
		}
		else
		{
			diff = (CDF-u) / CDF1;
		}	
        t = t - diff; 
        
        if (t < lower || t > upper)
        {
            t = 0.5 * (lower + upper);
        }

        CDF1 = dMixKT(t, v, par, pC, C1, C2, dC);
        CDF = C2MixKT(t, v, par, pC, C2);

        if ( CDF > u )
        {
            upper = t;
        }
        else
        {
            lower = t;
        }
    }
    return t;
}

double dKT(double u, double v, vector<double> par, Cfun pC, Cfun C1, Cfun C2, Cfun dC)
{
    // par:     parameters for the original copula, and a1 a2 for KT
    // pC:      cdf of the original copula
    // C1:      derivative to u1 of the original copula
    // C2:      derivative to u2 of the original copula
    // dC:      den of the original copula
 
    vector<double>::size_type nparC = par.size() - 2;     // nparC:   number of parameters of the original copula
    vector<double> parC(par.begin(), par.begin()+nparC );
    double a1 = par[nparC];
    double a2 = par[nparC+1];
    double ua1 = pow(u, a1);
    double va2 = pow(v, a2);
    double tem1 = (1-a1)*(1-a2) * pow(ua1, -1) * pow(va2, -1) * pC(ua1, va2, parC);
    double tem2 = (1-a1) * a2 * pow(ua1, -1) * C2(ua1, va2, parC);
    double tem3 = (1-a2) * a1 * pow(va2, -1) * C1(ua1, va2, parC);
    double tem4 = a1 * a2 * dC(ua1, va2, parC);
    double out = tem1 + tem2 + tem3 + tem4;
    return out;
}

double C1MixKT(double u, double v, vector<double> par, Cfun pC, Cfun C1)
{
    double out = 0.5 * (C1KT(u, v, par, pC, C1) - C1KT(1-u, 1-v, par, pC, C1) + 1);
    return out;
}

double C2MixKT(double u, double v, vector<double> par, Cfun pC, Cfun C2)
{
    double out = 0.5 * (C2KT(u, v, par, pC, C2) - C2KT(1-u, 1-v, par, pC, C2) + 1);
    return out;
}

// [[Rcpp::export]]
double C1MixKB5(double u, double v, vector<double> par)
{
    double out = C1MixKT(u, v, par, pB5, C1B5);
    return out;
}

// [[Rcpp::export]]
double C2MixKB5(double u, double v, vector<double> par)
{
    double out = C2MixKT(u, v, par, pB5, C2B5);
    return out;
}

#ifdef STUDENT_T
double C1MixKSt(double u, double v, vector<double> par)
{
    double out = C1MixKT(u, v, par, pSt, C1St);
    return out;
}

double C2MixKSt(double u, double v, vector<double> par)
{
    double out = C2MixKT(u, v, par, pSt, C2St);
    return out;
}
#endif

double dMixKT(double u, double v, vector<double> par, Cfun pC, Cfun C1, Cfun C2, Cfun dC)
{
    double out = 0.5 * (dKT(u, v, par, pC, C1, C2, dC) + dKT(1-u, 1-v, par, pC, C1, C2, dC));
    return out;
}

// [[Rcpp::export]]
double dMixKRB4(double u, double v, vector<double> par)
{
    double out = dMixKT(u, v, par, pRB4, C1RB4, C2RB4, dRB4);
    return out;
}

// [[Rcpp::export]]
double dMixKB4(double u, double v, vector<double> par)
{
    double out = dMixKT(u, v, par, pB4, C1B4, C2B4, dB4);
    return out;
}


// [[Rcpp::export]]
double dMixKB5(double u, double v, vector<double> par)
{
    double out = dMixKT(u, v, par, pB5, C1B5, C2B5, dB5);
    return out;
}

// [[Rcpp::export]]
double dMixKB6(double u, double v, vector<double> par)
{
    double out = dMixKT(u, v, par, pB6, C1B6, C2B6, dB6);
    return out;
}

#ifdef STUDENT_T
double dMixKSt(double u, double v, vector<double> par)
{
    double out = dMixKT(u, v, par, pSt, C1St, C2St, dSt);
    return out;
}
#endif

// based on qcondbb1 from CopulaModel R package
// [[Rcpp::export]]
double invC2B6(double u, double v, vector<double> par)
{ 
    double de,z,g,gp,x,y,con,de1,dif,mxdif,eps;
    int iter,mxiter;
    de = par[0];
    mxiter=20; 
    eps=1.e-6;
    x=-log(v); de1=de-1.;
    con=log(u)-x-de1*log(x); 
    z=x*pow(2.,1./de); // this is for z 
    mxdif=1; iter=0; 
    dif=.1;  // needed in case first step leads to NaN
    while(mxdif>eps && iter<mxiter)
    { 
	g=z+de1*log(z)+con;
	gp=1.+de1/z;
	if(isnan(g) || isnan(gp) || isnan(g/gp) )
        {
            dif/=-2.; 
	}  // added for de>50
	else dif=g/gp;
	z-=dif; 
	++iter;
	while(z<=x)
	{
		dif/=2.; z+=dif; 
	}
	mxdif=fabs(dif);
    }
    y=pow(pow(z,de)-pow(x,de),1./de);
    return exp(-y);
}

// [[Rcpp::export]]
double invC2BB1(double u, double v, vector<double> par)
{ 
    double de1,th1,vt,x,y,den,pden,diff,sm,smd,eps,r;
    double G21,gpdf,tem,thr;
    int iter,mxiter;
    double th,de;

    de=par[0]; 
    th=par[1];

    mxiter=30; 
    eps=1.e-6;
    de1=1./de; 
    th1=1./th;
    vt=pow(v,-th)-1.; 
    x=pow(vt,de); 
    den=pow(1.+vt,-th1-1)*vt/x; // density of univariate margin
    pden=u*den; // term 1/(th*de) cancels
    y=pow(pden,-1./(th1*de1+1))-x; // starting guess for y
    // in cases where x and ystart large pcond(qcond(q))-q might be .001
    if(x<1.e-5) // empirically based
    {
        y=x*(pow(u,-de/(de-1))-1.);  // take x=~ pow(th*(1-v),de) with v near 1?
        // need bound to prevent y large 
        if(y>1.e-5) y=1.e-5;
    }
	
    if(x>1.e5) // empirically based
    { 
        r=pow(u,-de*th/(1.+de*th))-1.;
        y=r*x;
        eps*=.0001; // good enough
            // some cases left where pcond(qcond(q))-q might be .001
    }
            // previous tolerance of 0.01
            // if th<0.1 or de<1.1 use boundary cases of Gumbel and MTCJ as starting points
            // *** biggest problem for de>5 and th small
            // v=(y^(1/de)+1)^(-1/th) so y=(v^(-th)-1)^de
    if(de<=1.1)   // MTCJ boundary
    { //qcondmtcj(&u,&v,&th,&v); y=pow(pow(v,-th)-1.,de); 
        thr=-th/(1.+th); tem=pow(u,thr)-1.;
        tem=tem*pow(v,-th)+1.; y=pow(tem-1.,de);
    }
	// Gumbel boundary if th<.2;
    else if(th<0.2)
    { 
        v = invC2B6(u,v,par);
        y = pow(pow(v,-th)-1., de);
    }

    diff=1; 
    iter=0;
    while(fabs(diff/y)> eps && iter<mxiter)
    { 
        sm=x+y; 
        smd=pow(sm,de1);
        G21=pow(1.+smd,-th1-1.) * smd/sm; 
        gpdf=-G21; 
        gpdf=gpdf/(1.+smd)/sm/de/th;
        gpdf=gpdf*(th*(de-1.)+(th*de+1.)*smd);
        iter=iter+1;
        diff=(G21-pden)/gpdf;
        y=y-diff;
        while(y<=0.)
        { 
            diff=diff/2.; 
            y=y+diff; 
        }
    }
    return pow(pow(y,de1)+1.,-th1);
}



// [[Rcpp::export]]
double invC2MixKB5(double u, double v, vector<double> par)
{
    double out = invC2MixKT(u, v, par, pB5, C1B5, C2B5, dB5);
    return out;
}

// [[Rcpp::export]]
double invC1MixKB5(double u, double v, vector<double> par)
{
    double out = invC1MixKT(u, v, par, pB5, C1B5, C2B5, dB5);
    return out;
}

#ifdef STUDENT_T
double invC1MixKSt(double u, double v, vector<double> par)
{
    double out = invC1MixKT(u, v, par, pSt, C1St, C2St, dSt);
    return out;
}
#endif


//// [[Rcpp::export]]
//double C1KB4test(double u, double v, vector<double> par)
//{
//    return C1KT(u, v, par, pB4, C1B4);
//}
//
//// [[Rcpp::export]]
//double C2KB4test(double u, double v, vector<double> par)
//{
//    return C2KT(u, v, par, pB4, C2B4);
//}


//////////////////////////////////////
// negative log likelihood of vines //
//////////////////////////////////////

// [[Rcpp::export]]
double nllk_vine(Rcpp::NumericVector par, Rcpp::NumericVector parM, Rcpp::NumericMatrix umat, const int Cs_family, const int C_family, const int mode, bool isPrint)
{
    // par: copula parameters
    // parM: parameter management, i.e. how many parameters of each
    // umat: matrix of the uniform scores data
    // Cs_family: copula families for C*
    // C_family: copula families for the exchangeable copulas
    // mode: 2: truncated  
    //       3: only 12 23 34 pairs
    //       4: only 13 24 pairs
    // isPrint: true for output details
    
    vector<double> cpar12, cpar23, cpar34, cpar13, cpar24;
    double exp_a1_12, exp_a2_12, exp_a1_23, exp_a2_23, exp_a1_34, exp_a2_34;
    int i;

    for(i = 0; i < parM[0]; ++i)
    {
        cpar12.push_back(par[i]);
    }

    for(i = 0; i < parM[1]; ++i)
    {
        cpar23.push_back(par[ i+parM[0] ]);
    }

    for(i = 0; i < parM[2]; ++i)
    {
        cpar34.push_back(par[ i+parM[0]+parM[1] ]);
    }

    for(i = 0; i < parM[3]; ++i)
    {
        cpar13.push_back(par[ i+parM[0]+parM[1]+parM[2] ]);
    }

    for(i = 0; i < parM[4]; ++i)
    {
        cpar24.push_back(par[ i+parM[0]+parM[1]+parM[2]+parM[3] ]);
    }

    if(Cs_family == 12 || Cs_family == 11 || Cs_family == -11) // 11: BB1,  12: Student t has two parameters (nu, r)
    {
        exp_a1_12 = exp(cpar12[2]);
        exp_a2_12 = exp(cpar12[3]);
        exp_a1_23 = exp(cpar23[2]);
        exp_a2_23 = exp(cpar23[3]);
        exp_a1_34 = exp(cpar34[2]);
        exp_a2_34 = exp(cpar34[3]);        
        
        // a1 and a2 should be between (0,1)
        cpar12[2] = exp_a1_12 / (exp_a1_12 + 1);
        cpar12[3] = exp_a2_12 / (exp_a2_12 + 1);
        cpar23[2] = exp_a1_23 / (exp_a1_23 + 1);
        cpar23[3] = exp_a2_23 / (exp_a2_23 + 1);
        cpar34[2] = exp_a1_34 / (exp_a1_34 + 1);
        cpar34[3] = exp_a2_34 / (exp_a2_34 + 1);
    }
    else
    {
        exp_a1_12 = exp(cpar12[1]);
        exp_a2_12 = exp(cpar12[2]);
        exp_a1_23 = exp(cpar23[1]);
        exp_a2_23 = exp(cpar23[2]);
        exp_a1_34 = exp(cpar34[1]);
        exp_a2_34 = exp(cpar34[2]);

        // a1 and a2 should be between (0,1)
        cpar12[1] = exp_a1_12 / (exp_a1_12 + 1);
        cpar12[2] = exp_a2_12 / (exp_a2_12 + 1);
        cpar23[1] = exp_a1_23 / (exp_a1_23 + 1);
        cpar23[2] = exp_a2_23 / (exp_a2_23 + 1);
        cpar34[1] = exp_a1_34 / (exp_a1_34 + 1);
        cpar34[2] = exp_a2_34 / (exp_a2_34 + 1);
    }
    
    int S = umat.nrow();  // sample size
    
    Rcpp::NumericVector uvec_1 = umat(Rcpp::_ , 0);
    Rcpp::NumericVector uvec_2 = umat(Rcpp::_ , 1);
    Rcpp::NumericVector uvec_3 = umat(Rcpp::_ , 2);
    Rcpp::NumericVector uvec_4 = umat(Rcpp::_ , 3);
    
	// Rcpp::Rcout << "dimension of data: " << S <<","<< umat.ncol() <<","<< uvec_1 <<","<< uvec_2 <<","<< uvec_3 <<","<< uvec_4 <<","<< std::endl;


    //-----------------------------------------------//
    //         initializing copula functions         //
    //-----------------------------------------------//

    // C1, C2, dC: for exchangeable copulas
    // Cs1, Cs2, pCs, dCs: for non-exchangeable copulas
    Cfun C2, dC, Cs1, Cs2, pCs, dCs;
    
    // Initializing parameters and choosing functions
    switch( Cs_family )
    {
        case -4:
        // codes for mixed KRB4 copula
            cpar12[0] = exp(cpar12[0]);
            cpar23[0] = exp(cpar23[0]);
            cpar34[0] = exp(cpar34[0]);
            pCs = &pRB4;
            dCs = &dRB4;
            Cs1 = &C1RB4;
            Cs2 = &C2RB4;
            break;
        case 4:
        // codes for mixed KB4 copula
            cpar12[0] = exp(cpar12[0]);
            cpar23[0] = exp(cpar23[0]);
            cpar34[0] = exp(cpar34[0]);
            pCs = &pB4;
            dCs = &dB4;
            Cs1 = &C1B4;
            Cs2 = &C2B4;
            break;
        case -5:
        // code for mixed KRB5 copula
            cpar12[0] = exp(cpar12[0]) + 1;
            cpar23[0] = exp(cpar23[0]) + 1;
            cpar34[0] = exp(cpar34[0]) + 1;
            pCs = &pRB5;
            dCs = &dRB5;
            Cs1 = &C1RB5;
            Cs2 = &C2RB5;
            break;
        case 5:
        // code for mixed KB5 copula
            cpar12[0] = exp(cpar12[0]) + 1;
            cpar23[0] = exp(cpar23[0]) + 1;
            cpar34[0] = exp(cpar34[0]) + 1;
            pCs = &pB5;
            dCs = &dB5;
            Cs1 = &C1B5;
            Cs2 = &C2B5;
            break;
        case -6:
        // code for mixed KRB6 copula
            cpar12[0] = exp(cpar12[0]) + 1;
            cpar23[0] = exp(cpar23[0]) + 1;
            cpar34[0] = exp(cpar34[0]) + 1;
            pCs = &pRB6;
            dCs = &dRB6;
            Cs1 = &C1RB6;
            Cs2 = &C2RB6;
            break;
        case 6:
        // code for mixed KB6 copula
            cpar12[0] = exp(cpar12[0]) + 1;
            cpar23[0] = exp(cpar23[0]) + 1;
            cpar34[0] = exp(cpar34[0]) + 1;
            pCs = &pB6;
            dCs = &dB6;
            Cs1 = &C1B6;
            Cs2 = &C2B6;
            break;
        case -11:
        // code for mixed KBB1 copula
            cpar12[0] = exp(cpar12[0]) + 1.;    // de >= 1, th > 0
            cpar12[1] = exp(cpar12[1]);
            cpar23[0] = exp(cpar23[0]) + 1.;
            cpar23[1] = exp(cpar23[1]);
            cpar34[0] = exp(cpar34[0]) + 1.;
            cpar34[1] = exp(cpar34[1]);
            pCs = &pRBB1;
            dCs = &dRBB1;
            Cs1 = &C1RBB1;
            Cs2 = &C2RBB1;
            break;
        case 11:
        // code for mixed KBB1 copula
            cpar12[0] = exp(cpar12[0]) + 1.;    // de >= 1, th > 0
            cpar12[1] = exp(cpar12[1]);
            cpar23[0] = exp(cpar23[0]) + 1.;
            cpar23[1] = exp(cpar23[1]);
            cpar34[0] = exp(cpar34[0]) + 1.;
            cpar34[1] = exp(cpar34[1]);
            pCs = &pBB1;
            dCs = &dBB1;
            Cs1 = &C1BB1;
            Cs2 = &C2BB1;
            break;
#ifdef STUDENT_T        
		case 12:
        // code for mixed KSt copula
            cpar12[0] = exp(cpar12[0]) + 1.;    // to-do: let nu > 1 to avoid extremal cases, maybe relaxed later
            cpar12[1] = (exp(2*cpar12[1]) - 1) / (exp(2*cpar12[1]) + 1);  // rho is between -1 and 1
            cpar23[0] = exp(cpar23[0]) + 1.;
            cpar23[1] = (exp(2*cpar23[1]) - 1) / (exp(2*cpar23[1]) + 1);
            cpar34[0] = exp(cpar34[0]) + 1.;
            cpar34[1] = (exp(2*cpar34[1]) - 1) / (exp(2*cpar34[1]) + 1);
            pCs = &pSt;
            dCs = &dSt;
            Cs1 = &C1St;
            Cs2 = &C2St;
            break;            
#endif			
        default:
          // code for mixed KB5 copula
            cpar12[0] = exp(cpar12[0]) + 1;
            cpar23[0] = exp(cpar23[0]) + 1;
            cpar34[0] = exp(cpar34[0]) + 1;
            pCs = &pB5;
            dCs = &dB5;
            Cs1 = &C1B5;
            Cs2 = &C2B5;
            break;
    }
    
    switch( C_family )
    {
        case 1:
        // codes for BB1 copula
            // BB1: de>=1, th>0
            cpar13[0] = exp(cpar13[0]) + 1;
            cpar13[1] = exp(cpar13[1]);
            cpar24[0] = exp(cpar24[0]) + 1;
            cpar24[1] = exp(cpar24[1]);
            dC = &dBB1;
            C2 = &C2BB1; // for exchangeable copula C1 and C2 are the same
            break;
        case -1:
        // codes for RT-BB1 copula
            // BB1: de>=1, th>0
            cpar13[0] = exp(cpar13[0]) + 1;
            cpar13[1] = exp(cpar13[1]);
            cpar24[0] = exp(cpar24[0]) + 1;
            cpar24[1] = exp(cpar24[1]);
            dC = &dRBB1;
            C2 = &C2RBB1; // for exchangeable copula C1 and C2 are the same
            break;    
        case 7:
        // codes for BB7 copula
            // BB7:  de>0, th>=1            
            cpar13[0] = exp(cpar13[0]);
            cpar13[1] = exp(cpar13[1]) + 1;
            cpar24[0] = exp(cpar24[0]);
            cpar24[1] = exp(cpar24[1]) + 1;
            dC = &dBB7;
            C2 = &C2BB7;
            break;
        default:
          // BB1: de>=1, th>0
            cpar13[0] = exp(cpar13[0]) + 1;
            cpar13[1] = exp(cpar13[1]);
            cpar24[0] = exp(cpar24[0]) + 1;
            cpar24[1] = exp(cpar24[1]);
            dC = &dBB1;
            C2 = &C2BB1; // for exchangeable copula C1 and C2 are the same
            break;
    }

    double u_1_2, u_3_2, u_2_3, u_4_3;
    double _ll12, _ll23, _ll34, _ll13_2, _ll24_3;
    
    double ll12 = 0;
    double ll23 = 0;
    double ll34 = 0;
    double ll13_2 = 0;
    double ll24_3 = 0;
		
    int Nfinite = 0; // number of records that have finite llk
    
    for(i = 0; i < S; ++i)
    {
        _ll12 = log(dMixKT((double)uvec_1[i], (double)uvec_2[i], cpar12, pCs, Cs1, Cs2, dCs));
        _ll23 = log(dMixKT((double)uvec_2[i], (double)uvec_3[i], cpar23, pCs, Cs1, Cs2, dCs));
        _ll34 = log(dMixKT((double)uvec_3[i], (double)uvec_4[i], cpar34, pCs, Cs1, Cs2, dCs));
        
        u_1_2 = C2MixKT((double)uvec_1[i], (double)uvec_2[i], cpar12, pCs, Cs2);
        u_3_2 = C1MixKT((double)uvec_2[i], (double)uvec_3[i], cpar23, pCs, Cs1);
        u_2_3 = C2MixKT((double)uvec_2[i], (double)uvec_3[i], cpar23, pCs, Cs2);
        u_4_3 = C1MixKT((double)uvec_3[i], (double)uvec_4[i], cpar34, pCs, Cs1);
        
        _ll13_2 = log(dC(u_1_2, u_3_2, cpar13));
        _ll24_3 = log(dC(u_2_3, u_4_3, cpar24));
                
        if(R_finite(_ll12) && R_finite(_ll23) && R_finite(_ll34) && R_finite(_ll13_2) && R_finite(_ll24_3))
        {
            ll12 += _ll12;
            ll23 += _ll23;
            ll34 += _ll34;
            ll13_2 += _ll13_2;
            ll24_3 += _ll24_3;
            ++Nfinite;
        }
        else
        {
            return 0;
        }
        
    }
    
    double out = 0;
    
    switch( mode )
    {
        case 2:
            out = -(ll12 + ll23 + ll34 + ll13_2 + ll24_3);
            break;
        case 3:
            out = -(ll12 + ll23 + ll34);
            break;
        case 4:
            out = -(ll13_2 + ll24_3);
            break;
        default:
            out = -(ll12 + ll23 + ll34 + ll13_2 + ll24_3);
            break;            
    }

    if(isPrint == true)
    {
        Rcpp::Rcout << "S/Nfinite/ll12/ll23/ll34/ll13_2/ll24_3: " << S << " " << Nfinite << " " << ll12 << " " << ll23 << " " << ll34 << " " << ll13_2 << " " << ll24_3 << " " << " " << out << " " << std::endl;
    }
    
    if(R_finite(out))
    {
        return out;
    }
    else
    {
        return 0;
    }
}

