#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

void d_f_j1_g(double a0, double a2, double b0, double b2, double q0, double q2, double r0, double tij, double t0, double geno_a, double geno_b, double geno_q, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;

	double e1 = pow(b, 2);
	double e4 = sqrt(2 * (e1 * q) + pow(a, 2));
	double e5 = 2 * q;
	double e6 = a + e4;
	double e7 = e6 / e5;
	double e8 = r0 - e7;
	double e9 = tij - t0;
	double e12 = 2 * (e4 / e8) + e5;
	double e13 = exp(2 * (e4 * e9));
	double e15 = e12 * e13 - e5;
	double e16 = a / e4;
	double e17 = e7 + 2 * (e4 / e15);
	double e20 = 1 + e16;
	double e21 = 2 * (e6 / pow(e5, 2));
	double e22 = 2 * e16;
	double e23 = 2 * (e1 / e4);
	double e25 = 4 * (q / e4);
	double e26 = e1 / (2 * (q * e4));

	// double df[3];
	double d_a = (-1)*(0.5 * ((e20 / e5 + (e22 - 2 * (((e20 * e4 / (q * e8) + e22) / e8 + 2 * (a * e12 * e9 / e4)) * e13 *e4 / e15)) / e15) / e17));
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = (-1)*(0.5 * (b * ((e25 - 2 * (((2 / e8 + e25) / e8 + 4 * (q * e12 * e9 / e4)) * e13 *e4 / e15)) / e15 + 1 / e4) / e17));
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = (-1)*(0.5 * (((e23 - 2 * ((((2 * ((e26 - e21) * e4 / e8) + e23) / e8 + 2 + 2 * (e1 * e12 * e9 / e4)) * e13 - 2) * e4 / e15)) / e15 + e26 - e21) / e17));
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;

	// return(df);
}

void d_f_j2_g(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f1, double m0, double r0, double tij, double yij, double t0, double geno_a, double geno_b, double geno_q, double geno_f, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;

	double e1 = pow(a, 2);
	double e2 = pow(b, 2);
	double e3 = e2 * q;
	double e5 = 2 * e3 + e1;
	double e6 = sqrt(e5);
	double e7 = 2 * q;
	double e8 = a + e6;
	double e9 = e8 / e7;
	double e10 = r0 - e9;
	double e11 = tij - t0;
	double e13 = 2 * (e6 / e10) + e7;
	double e14 = f - f1;
	double e16 = exp(2 * (e6 * e11));
	double e17 = e13 * e16;
	double e18 = a * e14;
	double e19 = e18 / e6;
	double e20 = e1 * e14;
	double e24 = 2 * (e20 / e6) - 2 * ((f - m0) * e6);
	double e26 = e7 - e17;
	double e27 = exp((-1)*t0 * e6);
	double e28 = exp(tij * e6);
	double e29 = e17 - e7;
	double e30 = e19 + f;
	double e39 = q * (2 * (e24 * e27 * e28 / e8) - 2 * (e19 - f)) - e13 * e30 * e16;
	double e40 = e39 / e26;
	double e41 = a / e6;
	double e42 = e9 + 2 * (e6 / e29);
	double e43 = yij - e40;
	double e46 = 1 + e41;
	double e47 = 2 * e8 / pow(e7, 2);
	double e49 = 2 * e6;
	double e50 = e2 / (2 * (q * e6));
	double e51 = pow(e49, 2);
	double e52 = 2 * (e2 / e6);
	double e53 = e50 - e47;
	double e56 = (2 * (e53 * e6 / e10) + e52) / e10 + 2 + 2 * (e2 * e13 * e11 / e6);
	double e57 = 0.5 * (t0 * e24);
	double e58 = 2 * e41;
	double e59 = 4 * (q / e6);
	double e60 = e1 / e5;
	double e62 = (e46 * e6 / (q * e10) + e58) / e10 + 2 * (a * e13 * e11 / e6);
	double e63 = e56 * e16;
	double e64 = e42 * e26;
	double e65 = e51 * e6;
	double e67 = (2 / e10 + e59) / e10 + 4 * (q * e13 * e11 / e6);
	double e68 = e40 - e30;
	double e71 = e57 + 4 * (e20 / e51) + f - m0;

	//double df[5];
	double d_a = ((e62 * e68 * e16 + q * (2 * (a * ((8 - 4 * e60) * e14 / 2 - 2 * (e57 + f - m0)) / e6 - e46 * e24 / e8) + 2 * (a * tij * e24 / e6)) * e27 * e28 / e8 - ((2 - 2 * e60) * e13 * e16 / 2 + 2 * (q * (1 - e60))) * e14 / e6) / e26 + 0.5 * ((e46 / e7 + (e58 - 2 * (e62 * e16 * e6 / e29)) / e29) * e43 / e42)) * e43 / e42;
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = b * (((e67 * e68 + 8 * (a * q * e13 * e14 / e65)) * e16 + pow(q, 2) * ((4 * (tij * e24) - 2 * (2 * (e24 / e8) + 4 * e71)) * e27 * e28 / e8 + 4 * (e18 / e5)) / e6) / e26 + 0.5 * (((e59 - 2 * (e67 * e16 * e6 / e29)) / e29 + 1 / e6) * e43 / e42)) * e43 / e42;
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = (((4 * (a * e2 * e13 * e14 / e65) - e56 * e30) * e16 + q * (2 * (e2 * tij * e24 / e6) - 2 * (2 * (e2 * e71 / e6) + 2 * (q * e24 * e53 / e8))) * e27 * e28 / e8 - ((2 - e63) * e39 / e26 + 2 * (a * (1 - e3 / e5) * e14 / e6 - f))) / e26 + 0.5 * (((e52 - 2 * ((e63 - 2) * e6 / e29)) / e29 + e50 - e47) * e43 / e42)) * e43 / e42;
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;
	double d_f = (q * (2 * ((2 * (e1 / e6) - e49) * e27 * e28 / e8) - 2 * (e41 - 1)) - e46 * e13 * e16) * e43 / e64;
	df[6] = d_f*(1 - 0.5*geno_f);
	df[7] = d_f*0.5*geno_f;
	df[8] = a * (e17 + q * (2 - 4 * (a * e27 * e28 / e8))) * e43 / (e64 * e6);

	//return(df);

}

void d_f_i1_g(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f1, double mu0, double theta, double m0, double r0, double tau, double t0, double geno_a, double geno_b, double geno_q, double geno_f, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;

	double e1 = pow(a, 2);
	double e2 = pow(b, 2);
	double e3 = e2 * q;
	double e5 = 2 * e3 + e1;
	double e6 = sqrt(e5);
	double e7 = 2 * q;
	double e8 = a + e6;
	double e10 = r0 - e8 / e7;
	double e11 = tau - t0;
	double e14 = 2 * (e6 / e10) + e7;
	double e15 = exp(2 * (e6 *e11));
	double e17 = e7 - e14 *e15;
	double e18 = f - f1;
	double e19 = e1 *e18;
	double e23 = 2 * (e19 / e6) - 2 * ((f - m0) *e6);
	double e25 = exp((-1)*t0 *e6);
	double e26 = exp(tau *e6);
	double e27 = e17 *e8;
	double e28 = 1 - 4 * (q / e17);
	double e29 = a *e28;
	double e31 = e29 *e18 / e6;
	double e32 = q *e23;
	double e34 = 2 * (e32 *e25 *e26 / e27) + e31;
	double e35 = exp(tau * theta);
	double e36 = a / e6;
	double e41 = 0.5 *e8 + mu0 *e35 + q * (pow(e34, 2) - 2 * (e6 / e17));
	double e42 = 2 * e6;
	double e43 = 1 + e36;
	double e49 = e2 / (2 * (q *e6)) - 2 * (e8 / pow(e7, 2));
	double e50 = e2 / e6;
	double e52 = pow((-1)*e17, 2);
	double e54 = (e43 *e6 / (q *e10) + 2 * e36) / e10 + 2 * (a *e14 *e11 / e6);
	double e56 = ((2 / e10 + 4 * (q / e6)) / e10 + 4 * (q *e14 *e11 / e6)) *e15;
	double e57 = 0.5 * t0 *e23;
	double e58 = 2 - ((2 * (e49 *e6 / e10) + 2 * e50) / e10 + 2 + 2 * (e2 *e14 *e11 / e6)) *e15;
	double e59 = e54 *e15;
	double e60 = e17 *e6;
	double e61 = pow(e17, 2);
	double e65 = e57 + 4 * (e19 / pow(e42, 2)) + f - m0;
	double e66 = a * q;
	double e67 = e1 / e5;
	double e68 = q *e58;

	// double df[11];
	double d_a = (0.5 *e43 + q * (2 * (((e28 * (2 - 2 * e67) - 8 * (e66 *e54 *e15 / e61)) *e18 / e42 + q * (2 * ((e59 / e17 + a * tau / e6) *e23) + 2 * (a * ((8 - 4 * e67) *e18 / 2 - 2 * (e57 + f - m0)) / e6 - e43 *e23 / e8)) *e25 *e26 / e27) *e34) - (2 * (e59 *e6 / e52) + 2 * (a / e60)))) / e41;
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = b * q * (1 / e6 + q * (2 * (((2 * ((e56 / e17 + 2 * (q * tau / e6)) *e23) - 2 * (q * (2 * (e23 / e8) + 4 * e65) / e6)) *e25 *e26 / e27 - a * (4 * (e28 / e5) + 8 * (e56 / e61)) *e18 / e42) *e34) - 4 / e60) - 2 * (e56 *e6 / e52)) / e41;
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = (e34 * (e31 + q * (2 * (e23 *e25 *e26 / e27) + 2 * (q * (2 * (e23 * (e2 * tau / e6 - e58 / e17)) - 2 * (2 * (e2 *e65 / e6) + 2 * (e32 *e49 / e8))) *e25 *e26 / e27 - a * (2 * ((4 - 4 * (e68 / e17)) / e17) + 2 * (e2 *e28 / e5)) *e18 / e42))) + 0.5 *e50 + 2 * (e68 *e6 / e52) - (2 * (e3 / e6) + e42) / e17) / e41;
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;
	double d_f = 2 * (q *e34 * (2 * (q * (2 * (e1 / e6) - e42) *e25 *e26 / e27) + e29 / e6) / e41);
	df[6] = d_f*(1 - 0.5*geno_f);
	df[7] = d_f*0.5*geno_f;
	df[8] = -(2 * (e66 *(1 + q * (4 * (a *e25 *e26 / e8) - 4) / e17) *e34 / (e41 *e6)));
	df[9] = e35 / e41;
	df[10] = mu0 * tau *e35 / e41;

	// return(df);
}

void dev_mu_int_g(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f1, double mu0, double theta, double tij, double r0, double m0, double t0, double geno_a, double geno_b, double geno_q, double geno_f, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;

	double e1 = pow(b, 2);
	double e2 = pow(a, 2);
	double e3 = e1 * q;
	double e5 = 2 * e3 + e2;
	double e6 = sqrt(e5);
	double e7 = 2 * q;
	double e8 = a + e6;
	double e10 = r0 - e8 / e7;
	double e11 = 2 * (e6 / e10);
	double e12 = e11 + e7;
	double e13 = t0 *e6;
	double e15 = exp(-2 * e13);
	double e16 = tij - t0;
	double e17 = e12 *e15;
	double e19 = sqrt(-0.5 * e17 / q);
	double e20 = 2 * (e6 *e16);
	double e21 = f - f1;
	double e22 = exp(e20);
	double e23 = 2 * e6;
	double e25 = e7 - e12 *e22;
	double e26 = e2 *e21;
	double e27 = e26 / e6;
	double e31 = 2 * e27 - 2 * ((f - m0) *e6);
	double e32 = a / e6;
	double e33 = tij *e6;
	double e34 = 2 * (q *e6);
	double e35 = pow(e23, 2);
	double e36 = pow(e7, 2);
	double e37 = exp(e13);
	double e38 = exp(e33);
	double e39 = 1 + e32;
	double e43 = e1 / e34 - 2 * (e8 / e36);
	double e44 = e1 / e6;
	double e48 = e39 *e6 / (q *e10) + 2 * e32;
	double e52 = 2 * (e43 *e6 / e10) + 2 * e44;
	double e54 = 2 / e10 + 4 * (q / e6);
	double e55 = exp((-1)*e13);
	double e56 = e37 *e19;
	double e57 = e38 *e19;
	double e59 = e12 *e35 *e15;
	double e60 = (-1)*e11;
	double e61 = e10 / e23;
	double e62 = 1 / e25;
	double e63 = e48 / e10;
	double e65 = e52 / e10 + 2;
	double e66 = e54 / e10;
	double e69 = atan(e57) - atan(e56);
	double e70 = log(e60);
	double e71 = log(e25);
	double e72 = q *e19;
	double e73 = a * t0;
	double e74 = e1 * t0;
	double e75 = 1 / e7;
	double e76 = q * t0;
	double e77 = e61 + e62;
	double e79 = e63 - 2 * (e73 *e12 / e6);
	double e80 = e65 - 2 * (e74 *e12 / e6);
	double e81 = e17 *e6;
	double e82 = e66 - 4 * (e76 *e12 / e6);
	double e83 = 0.5 * (t0 *e31);
	double e84 = a *e12;
	double e85 = e1 *e12;
	double e86 = exp((-1)*e33);
	double e87 = q *e12;
	double e88 = e71 - e70;
	double e89 = a *e21;
	double e90 = (e63 + 2 * (e84 *e16 / e6)) *e22;
	double e91 = (e65 + 2 * (e85 *e16 / e6)) *e22;
	double e92 = pow(e59, 2);
	double e93 = (e66 + 4 * (e87 *e16 / e6)) *e22;
	double e96 = (e61 + e75) *e55 + (e62 - e75) *e86 - e69 *e19 / e7;
	double e99 = pow(e56, 2) + 1;
	double e101 = pow(e57, 2) + 1;
	double e104 = 0.5 * (e12 / q) - 0.5 *e80;
	double e107 = e83 + 4 * (e26 / e35) + f - m0;
	double e109 = e2 / e5;
	double e110 = e71 - (e20 + e70);
	double e111 = pow(q, 2);
	double e113 = e77 / e23 + (e16 / 2 - e88 / (4 * e6)) / q;
	double e114 = e104 *e15;
	double e115 = pow(e89 / e23, 2);
	double e116 = pow(e8, 2);
	double e117 = pow(e55, 2);
	double e123 = pow(e23, 3);
	double e124 = 0.25 * (e79 *e15 / e72);
	double e125 = 0.25 * (e82 *e15 / e72);
	double e127 = 2 * e81;
	double e128 = 32 * e31;
	double e129 = a * tij;
	double e130 = e1 * tij;
	double e131 = q *e31;
	double e132 = q *e25;
	double e135 = pow(e60, 2) * e10;
	double e137 = pow((-4) * e5 / e10, 2) * e10;
	double e139 = e90 + e48 *e25 / e23;
	double e143 = e93 + e25 *e54 / e23;
	double e146 = pow(e127, 2);
	double e147 = pow(2 * e25 *e6, 2);
	double e149 = pow(e34, 2);
	double e151 = e12 *e8 *e15;
	double e152 = pow(e25, 2);
	double e156 = pow(e21, 2);
	double e157 = 16 * e31;
	double e158 = 2 - (e91 + e52 *e25 / e23);
	double e159 = 2 - e91;
	double e160 = 2 * (e79 *e6);
	double e161 = 2 * (e80 *e6);
	double e162 = 2 * (e82 *e6);
	double e164 = 2 * (e31 / e8) + 4 * e107;
	double e165 = 2 * e5;
	double e167 = 2 * (e2 / e6) - e23;
	double e169 = 2 * (e1 *e107 / e6) + 2 * (e131 *e43 / e8);
	double e172 = a * ((8 - 4 * e109) *e21 / 2 - 2 * (e83 + f - m0)) / e6 - e39 *e31 / e8;
	double e173 = e84 / e6;
	double e174 = a * q;
	double e175 = e85 / e6;
	double e176 = exp(t0 * theta);
	double e177 = exp(theta * tij);
	double e178 = q *e113;
	double e179 = e87 / e6;
	double e180 = q * tij;
	double e182 = e139 / e25 + 2 * (a *e16 / e6);
	double e184 = e143 / e25 + 4 * (q *e16 / e6);
	double e185 = ((0.5 * (e114 / e72) + e130 *e19 / e6) *e38 / e101 - (0.5 * (e114 / e72) + e74 *e19 / e6) *e37 / e99) *e19;
	double e187 = e158 / e25 - 2 * (e1 *e16 / e6);
	double e188 = e167 *e21;
	double e189 = e5 *e12;
	double e194 = (2 * (e180 *e19 / e6) - e125) *e38 / e101 - (2 * (e76 *e19 / e6) - e125) *e37 / e99;
	double e195 = e123 *e6;
	double e200 = (e129 *e19 / e6 - e124) *e38 / e101 - (e73 *e19 / e6 - e124) *e37 / e99;
	double e201 = (e177 - e176) / theta;
	double e202 = e16 / e6;
	double e203 = 1 - e109;
	double e204 = e160 + 4 * e173;
	double e206 = e162 + 8 * e179;
	double e207 = 2 / e36;
	double e208 = 4 * e132;
	double e209 = 64 * e178;
	double e210 = 8 * e72;
	double e211 = 8 * e16;
	double e212 = a *e172;
	double e214 = e89 *e110 / e6;
	double e215 = q *e96;
	double e216 = q * (e161 + 4 * e175);
	double e217 = q *e164;
	double e218 = q *e169;
	double e220 = e111 *e77 *e31;

	// double df[11];
	double d_a = 0.5 * (e39 *e16) + q * (((((4 * (a *e79 *e31 *e15 / e19) - q * (e157 + 16 * e212) *e19) / e59 + 32 * (e174 *e204 *e31 *e15 *e19 *e6 / e92)) *e69 - q * (16 * (a *e200 *e31 *e19 / e59) + q *(e96 * ((e128 + 32 * e212) / e59 - 64 * (a *e204 *e31 *e15 *e6 / e92)) + 32 * (a * ((e129 / e34 - (e129 / e6 - e90 / e25) / e25) *e86 - (((e48 / 2 + e73) *e10 / e165 + e73 / e34) *e55 + (e200 *e19 / 2 - e79 *e69 *e15 / e210) / q)) *e31 / e59)))) *e55 / e8 + a * (8 * (e203 *e16) - ((48 * e109 - 32) *e110 + 16 * (a *e182)) / e23) *e21 / e35) *e21 + q * (2 * (q * ((e90 / e152 - e48 / e135) *e31 / e81 + e77 * (2 * (e172 / e81) - 2 * ((e160 + 2 * e173) *e31 *e15 / e146))) *e31 *e117 / e116) + 32 * (((e139 / e208 + e174 *e88 / e149) / e6 - ((2 * (e48 *e6) + 4 * a) / e137 + (2 * (a *e25 / e6) - 2 * (e90 *e6)) / e147)) *e115) + 64 * (a *e113 *e203 *e156 / e35))) - 0.5 *e182;
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = b *(q * (e202 + a * ((((16 * (e111 *e164 *e19 / e6) + 4 * (e82 *e31 *e15 / e19)) / e59 + 32 * (q *e206 *e31 *e15 *e19 *e6 / e92)) *e69 - q * (16 * (e194 *e31 *e19 / e59) + q * (32 * (e31 * (e86 *(tij / e6 - (2 * (e180 / e6) - e93 / e25) / e25) - ((e194 *e19 / 2 - e82 *e69 *e15 / e210) / q + ((e54 / 2 + 2 * e76) *e10 / e165 + t0 / e6) *e55)) / e59) - e96 * (32 * (e217 / (e59 *e6)) + 64 * (e206 *e31 *e15 *e6 / e92))))) *e55 / e8 - a * (16 * e184 + 96 * (q *e110 / e5)) *e21 / e123) *e21 + q * (2 * (q * ((e93 / e152 - e54 / e135) *e31 / e81 - e77 * (2 * ((e162 + 4 * e179) *e31 *e15 / e146) + 2 * (e217 / (e189 *e15)))) *e31 *e117 / e116) + 32 * (((e143 / e208 + 2 * (e111 *e88 / e149)) / e6 - ((2 * (e54 *e6) + 8 * q) / e137 + (4 * (e132 / e6) - 2 * (e93 *e6)) / e147)) *e115) - e2 * (256 * e178 + 32 * e16) *e156 / e195)) - 0.5 *e184);
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = e77 * pow(2 * e131 *e55 / e8, 2) / e127 + (0.5 *e44 + 4 * e115) *e16 + 0.5 *e187 + q * (e113 *(64 * e115 - 128 * (e2 *e1 * q *e156 / e195)) + a * (((e31 * (32 * (e216 *e15 *e6 / e92) - 16 / e59) *e19 - (8 * e104 *e31 *e15 / e19 - 16 * (e218 *e19)) / e59) *e69 - q * (e96 * ((e128 - 32 * e218) / e59 - 64 * (e216 *e31 *e15 *e6 / e92)) + (16 * e185 + 32 * (q * ((e207 + e130 / e34 - (e159 / e25 + e130 / e6) / e25) *e86 - (((e52 / 2 + e74) *e10 / e165 + e207 + e74 / e34) *e55 + e185 / e7 + (e114 / (4 * e111 *e19) - 2 * (e19 / e36)) *e69)))) *e31 / e59)) *e55 / e8 - a * (e1 *(16 * e202 + 48 * (e110 / e5)) - 16 * e187) *e21 / e123) *e21 + q * (32 * ((0.5 * ((2 * (e3 / e6) + e23) *e88 / e149) - (e158 / (4 * (e132 *e6)) + (2 * (e159 *e6) + 2 * (e1 *e25 / e6)) / e147 + (2 * (e52 *e6) + 4 * e1) / e137 + 2 * (e16 / e36))) *e115) - 2 * (q *((e159 / e152 + e52 / e135) *e31 / e81 + e77 * (2 * ((e161 + 2 * e175) *e31 *e15 / e146) + 2 * (e169 / e81))) *e31 *e117 / e116))) - a * (32 * (e111 *e96 *e31 *e55 / e151) - 8 * e214) *e21 / e35;
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;
	double d_f = q * (4 * (e220 *e167 *e117 / (e12 *e116 *e15 *e6)) + a * (a *(16 * (e110 / e6) + e209 + e211) *e21 - q * ((16 * e188 + e157) *e69 *e19 + e215 * (32 * e188 + e128)) *e55 / e151) / e35);
	df[6] = d_f*(1 - 0.5*geno_f);
	df[7] = d_f*0.5*geno_f;
	df[8] = e174 * ((q * ((e157 + 32 * e27) *e69 *e19 + e215 * (e128 + 64 * e27)) *e55 / e151 - 16 * e214) / e35 - a * ((e209 + e211) *e21 / e35 + 8 * (e220 *e117 / (e189 *e116 *e15))));
	df[9] = e201;
	df[10] = 0;
	if (theta != 0)
	{
		df[10] = mu0 * (tij *e177 - (e201 + t0 *e176)) / theta;
	}


	// return(df);
}

RcppExport SEXP devlik_g(SEXP param, SEXP m0, SEXP r0, SEXP tau, SEXP yij, SEXP delta, SEXP tij, SEXP n_j, SEXP geno_a, SEXP geno_b, SEXP geno_q, SEXP geno_f, SEXP t0)
{

	arma::vec param_c = as<arma::vec>(param);
	// double lik = 0;

	arma::vec m0_v = as<arma::vec>(m0);
	arma::vec r0_v = as<arma::vec>(r0);
	arma::vec tau_v = as<arma::vec>(tau);
	arma::vec y_v = as<arma::vec>(yij);
	arma::vec delta_v = as<arma::vec>(delta);
	arma::vec t_v = as<arma::vec>(tij);
	arma::vec n_v = as<arma::vec>(n_j);
	arma::vec t0_v = as<arma::vec>(t0);
	arma::vec geno_a_v = as<arma::vec>(geno_a);
	arma::vec geno_b_v = as<arma::vec>(geno_b);
	arma::vec geno_q_v = as<arma::vec>(geno_q);
	arma::vec geno_f_v = as<arma::vec>(geno_f);

	double a0 = param_c(0);
	double a2 = param_c(1);
	double b0 = param_c(2);
	double b2 = param_c(3);
	double q0 = param_c(4);
	double q2 = param_c(5);
	double f0 = param_c(6);
	double f2 = param_c(7);
	double f1 = param_c(8);
	double mu0 = param_c(9);
	double theta = param_c(10);

	int i_n = n_v.n_elem;
	int t_n = m0_v.n_elem;
	int index = 0;

	double dev_a0 = 0;
	double dev_a2 = 0;
	double dev_b0 = 0;
	double dev_b2 = 0;
	double dev_q0 = 0;
	double dev_q2 = 0;
	double dev_f0 = 0;
	double dev_f2 = 0;
	double dev_f1 = 0;
	double dev_mu0 = 0;
	double dev_theta = 0;

	for (int i = 0; i < i_n; i++)
	{
		int start = index;
		int end = index + n_v(i);
		double geno_a_i = geno_a_v(i);
		double geno_b_i = geno_b_v(i);
		double geno_q_i = geno_q_v(i);
		double geno_f_i = geno_f_v(i);

		for (int j = start; j < end; j++)
		{
			double re_1[6];
			for (int k = 0; k < 6; k++)
			{
				re_1[k] = 0;
			}
			d_f_j1_g(a0, a2, b0, b2, q0, q2, r0_v[j], t_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, re_1);
			double re_2[9];
			for (int k = 0; k < 9; k++)
			{
				re_2[k] = 0;
			}
			d_f_j2_g(a0, a2, b0, b2, q0, q2, f0, f2, f1, m0_v[j], r0_v[j], t_v[j], y_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, geno_f_i, re_2);
			double re_4[11];
			for (int k = 0; k < 11; k++)
			{
				re_4[k] = 0;
			}

			dev_mu_int_g(a0, a2, b0, b2, q0, q2, f0, f2, f1, mu0, theta, t_v[j], r0_v[j], m0_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, geno_f_i, re_4);

			dev_a0 = dev_a0 + re_1[0] + re_2[0] - re_4[0];
			dev_a2 = dev_a2 + re_1[1] + re_2[1] - re_4[1];
			dev_b0 = dev_b0 + re_1[2] + re_2[2] - re_4[2];
			dev_b2 = dev_b2 + re_1[3] + re_2[3] - re_4[3];
			dev_q0 = dev_q0 + re_1[4] + re_2[4] - re_4[4];
			dev_q2 = dev_q2 + re_1[5] + re_2[5] - re_4[5];
			dev_f0 = dev_f0 + re_2[6] - re_4[6];
			dev_f2 = dev_f2 + re_2[7] - re_4[7];
			dev_f1 = dev_f1 + re_2[8] - re_4[8];
			dev_mu0 -= re_4[9];
			dev_theta -= re_4[10];
		}

		int pre_end = end - 1;
		double re_3[11];

		if (delta_v(i)>0)
		{
			d_f_i1_g(a0, a2, b0, b2, q0, q2, f0, f2, f1, mu0, theta, m0_v(pre_end), r0_v(pre_end), tau_v(i), t0_v(pre_end), geno_a_i, geno_b_i, geno_q_i, geno_f_i, re_3);

			dev_a0 += re_3[0];
			dev_a2 += re_3[1];
			dev_b0 += re_3[2];
			dev_b2 += re_3[3];
			dev_q0 += re_3[4];
			dev_q2 += re_3[5];
			dev_f0  += re_3[6];
			dev_f2 += re_3[7];
			dev_f1 += re_3[8];
			dev_mu0 += re_3[9];
			dev_theta += re_3[10];

		}

		index = end;
	}

	arma::vec res = arma::zeros<arma::vec>(11);
	res(0) = dev_a0*(-1);
	res(1) = dev_a2*(-1);
	res(2) = dev_b0*(-1);
	res(3) = dev_b2*(-1);
	res(4) = dev_q0*(-1);
	res(5) = dev_q2*(-1);
	res(6) = dev_f0*(-1);
	res(7) = dev_f2*(-1);
	res(8) = dev_f1*(-1);
	res(9) = dev_mu0*(-1);
	res(10) = dev_theta*(-1);

	return(Rcpp::wrap(res));

}
