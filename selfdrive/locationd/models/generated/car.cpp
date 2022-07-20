#include "car.h"

namespace {
#define DIM 9
#define EDIM 9
#define MEDIM 9
typedef void (*Hfun)(double *, double *, double *);

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}
const static double MAHA_THRESH_25 = 3.8414588206941227;
const static double MAHA_THRESH_24 = 5.991464547107981;
const static double MAHA_THRESH_30 = 3.8414588206941227;
const static double MAHA_THRESH_26 = 3.8414588206941227;
const static double MAHA_THRESH_27 = 3.8414588206941227;
const static double MAHA_THRESH_29 = 3.8414588206941227;
const static double MAHA_THRESH_28 = 3.8414588206941227;
const static double MAHA_THRESH_31 = 3.8414588206941227;

/******************************************************************************
 *                       Code generated with sympy 1.9                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_5872706819683017341) {
   out_5872706819683017341[0] = delta_x[0] + nom_x[0];
   out_5872706819683017341[1] = delta_x[1] + nom_x[1];
   out_5872706819683017341[2] = delta_x[2] + nom_x[2];
   out_5872706819683017341[3] = delta_x[3] + nom_x[3];
   out_5872706819683017341[4] = delta_x[4] + nom_x[4];
   out_5872706819683017341[5] = delta_x[5] + nom_x[5];
   out_5872706819683017341[6] = delta_x[6] + nom_x[6];
   out_5872706819683017341[7] = delta_x[7] + nom_x[7];
   out_5872706819683017341[8] = delta_x[8] + nom_x[8];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4273870982427227725) {
   out_4273870982427227725[0] = -nom_x[0] + true_x[0];
   out_4273870982427227725[1] = -nom_x[1] + true_x[1];
   out_4273870982427227725[2] = -nom_x[2] + true_x[2];
   out_4273870982427227725[3] = -nom_x[3] + true_x[3];
   out_4273870982427227725[4] = -nom_x[4] + true_x[4];
   out_4273870982427227725[5] = -nom_x[5] + true_x[5];
   out_4273870982427227725[6] = -nom_x[6] + true_x[6];
   out_4273870982427227725[7] = -nom_x[7] + true_x[7];
   out_4273870982427227725[8] = -nom_x[8] + true_x[8];
}
void H_mod_fun(double *state, double *out_1982852844921793486) {
   out_1982852844921793486[0] = 1.0;
   out_1982852844921793486[1] = 0;
   out_1982852844921793486[2] = 0;
   out_1982852844921793486[3] = 0;
   out_1982852844921793486[4] = 0;
   out_1982852844921793486[5] = 0;
   out_1982852844921793486[6] = 0;
   out_1982852844921793486[7] = 0;
   out_1982852844921793486[8] = 0;
   out_1982852844921793486[9] = 0;
   out_1982852844921793486[10] = 1.0;
   out_1982852844921793486[11] = 0;
   out_1982852844921793486[12] = 0;
   out_1982852844921793486[13] = 0;
   out_1982852844921793486[14] = 0;
   out_1982852844921793486[15] = 0;
   out_1982852844921793486[16] = 0;
   out_1982852844921793486[17] = 0;
   out_1982852844921793486[18] = 0;
   out_1982852844921793486[19] = 0;
   out_1982852844921793486[20] = 1.0;
   out_1982852844921793486[21] = 0;
   out_1982852844921793486[22] = 0;
   out_1982852844921793486[23] = 0;
   out_1982852844921793486[24] = 0;
   out_1982852844921793486[25] = 0;
   out_1982852844921793486[26] = 0;
   out_1982852844921793486[27] = 0;
   out_1982852844921793486[28] = 0;
   out_1982852844921793486[29] = 0;
   out_1982852844921793486[30] = 1.0;
   out_1982852844921793486[31] = 0;
   out_1982852844921793486[32] = 0;
   out_1982852844921793486[33] = 0;
   out_1982852844921793486[34] = 0;
   out_1982852844921793486[35] = 0;
   out_1982852844921793486[36] = 0;
   out_1982852844921793486[37] = 0;
   out_1982852844921793486[38] = 0;
   out_1982852844921793486[39] = 0;
   out_1982852844921793486[40] = 1.0;
   out_1982852844921793486[41] = 0;
   out_1982852844921793486[42] = 0;
   out_1982852844921793486[43] = 0;
   out_1982852844921793486[44] = 0;
   out_1982852844921793486[45] = 0;
   out_1982852844921793486[46] = 0;
   out_1982852844921793486[47] = 0;
   out_1982852844921793486[48] = 0;
   out_1982852844921793486[49] = 0;
   out_1982852844921793486[50] = 1.0;
   out_1982852844921793486[51] = 0;
   out_1982852844921793486[52] = 0;
   out_1982852844921793486[53] = 0;
   out_1982852844921793486[54] = 0;
   out_1982852844921793486[55] = 0;
   out_1982852844921793486[56] = 0;
   out_1982852844921793486[57] = 0;
   out_1982852844921793486[58] = 0;
   out_1982852844921793486[59] = 0;
   out_1982852844921793486[60] = 1.0;
   out_1982852844921793486[61] = 0;
   out_1982852844921793486[62] = 0;
   out_1982852844921793486[63] = 0;
   out_1982852844921793486[64] = 0;
   out_1982852844921793486[65] = 0;
   out_1982852844921793486[66] = 0;
   out_1982852844921793486[67] = 0;
   out_1982852844921793486[68] = 0;
   out_1982852844921793486[69] = 0;
   out_1982852844921793486[70] = 1.0;
   out_1982852844921793486[71] = 0;
   out_1982852844921793486[72] = 0;
   out_1982852844921793486[73] = 0;
   out_1982852844921793486[74] = 0;
   out_1982852844921793486[75] = 0;
   out_1982852844921793486[76] = 0;
   out_1982852844921793486[77] = 0;
   out_1982852844921793486[78] = 0;
   out_1982852844921793486[79] = 0;
   out_1982852844921793486[80] = 1.0;
}
void f_fun(double *state, double dt, double *out_8224871313240091039) {
   out_8224871313240091039[0] = state[0];
   out_8224871313240091039[1] = state[1];
   out_8224871313240091039[2] = state[2];
   out_8224871313240091039[3] = state[3];
   out_8224871313240091039[4] = state[4];
   out_8224871313240091039[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] - 9.8000000000000007*state[8] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_8224871313240091039[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_8224871313240091039[7] = state[7];
   out_8224871313240091039[8] = state[8];
}
void F_fun(double *state, double dt, double *out_5540635882771833240) {
   out_5540635882771833240[0] = 1;
   out_5540635882771833240[1] = 0;
   out_5540635882771833240[2] = 0;
   out_5540635882771833240[3] = 0;
   out_5540635882771833240[4] = 0;
   out_5540635882771833240[5] = 0;
   out_5540635882771833240[6] = 0;
   out_5540635882771833240[7] = 0;
   out_5540635882771833240[8] = 0;
   out_5540635882771833240[9] = 0;
   out_5540635882771833240[10] = 1;
   out_5540635882771833240[11] = 0;
   out_5540635882771833240[12] = 0;
   out_5540635882771833240[13] = 0;
   out_5540635882771833240[14] = 0;
   out_5540635882771833240[15] = 0;
   out_5540635882771833240[16] = 0;
   out_5540635882771833240[17] = 0;
   out_5540635882771833240[18] = 0;
   out_5540635882771833240[19] = 0;
   out_5540635882771833240[20] = 1;
   out_5540635882771833240[21] = 0;
   out_5540635882771833240[22] = 0;
   out_5540635882771833240[23] = 0;
   out_5540635882771833240[24] = 0;
   out_5540635882771833240[25] = 0;
   out_5540635882771833240[26] = 0;
   out_5540635882771833240[27] = 0;
   out_5540635882771833240[28] = 0;
   out_5540635882771833240[29] = 0;
   out_5540635882771833240[30] = 1;
   out_5540635882771833240[31] = 0;
   out_5540635882771833240[32] = 0;
   out_5540635882771833240[33] = 0;
   out_5540635882771833240[34] = 0;
   out_5540635882771833240[35] = 0;
   out_5540635882771833240[36] = 0;
   out_5540635882771833240[37] = 0;
   out_5540635882771833240[38] = 0;
   out_5540635882771833240[39] = 0;
   out_5540635882771833240[40] = 1;
   out_5540635882771833240[41] = 0;
   out_5540635882771833240[42] = 0;
   out_5540635882771833240[43] = 0;
   out_5540635882771833240[44] = 0;
   out_5540635882771833240[45] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_5540635882771833240[46] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_5540635882771833240[47] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5540635882771833240[48] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_5540635882771833240[49] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_5540635882771833240[50] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_5540635882771833240[51] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_5540635882771833240[52] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_5540635882771833240[53] = -9.8000000000000007*dt;
   out_5540635882771833240[54] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_5540635882771833240[55] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_5540635882771833240[56] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5540635882771833240[57] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5540635882771833240[58] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_5540635882771833240[59] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_5540635882771833240[60] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_5540635882771833240[61] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_5540635882771833240[62] = 0;
   out_5540635882771833240[63] = 0;
   out_5540635882771833240[64] = 0;
   out_5540635882771833240[65] = 0;
   out_5540635882771833240[66] = 0;
   out_5540635882771833240[67] = 0;
   out_5540635882771833240[68] = 0;
   out_5540635882771833240[69] = 0;
   out_5540635882771833240[70] = 1;
   out_5540635882771833240[71] = 0;
   out_5540635882771833240[72] = 0;
   out_5540635882771833240[73] = 0;
   out_5540635882771833240[74] = 0;
   out_5540635882771833240[75] = 0;
   out_5540635882771833240[76] = 0;
   out_5540635882771833240[77] = 0;
   out_5540635882771833240[78] = 0;
   out_5540635882771833240[79] = 0;
   out_5540635882771833240[80] = 1;
}
void h_25(double *state, double *unused, double *out_3133515492953738915) {
   out_3133515492953738915[0] = state[6];
}
void H_25(double *state, double *unused, double *out_6812056836489383815) {
   out_6812056836489383815[0] = 0;
   out_6812056836489383815[1] = 0;
   out_6812056836489383815[2] = 0;
   out_6812056836489383815[3] = 0;
   out_6812056836489383815[4] = 0;
   out_6812056836489383815[5] = 0;
   out_6812056836489383815[6] = 1;
   out_6812056836489383815[7] = 0;
   out_6812056836489383815[8] = 0;
}
void h_24(double *state, double *unused, double *out_9067730718024463655) {
   out_9067730718024463655[0] = state[4];
   out_9067730718024463655[1] = state[5];
}
void H_24(double *state, double *unused, double *out_4216295603543240908) {
   out_4216295603543240908[0] = 0;
   out_4216295603543240908[1] = 0;
   out_4216295603543240908[2] = 0;
   out_4216295603543240908[3] = 0;
   out_4216295603543240908[4] = 1;
   out_4216295603543240908[5] = 0;
   out_4216295603543240908[6] = 0;
   out_4216295603543240908[7] = 0;
   out_4216295603543240908[8] = 0;
   out_4216295603543240908[9] = 0;
   out_4216295603543240908[10] = 0;
   out_4216295603543240908[11] = 0;
   out_4216295603543240908[12] = 0;
   out_4216295603543240908[13] = 0;
   out_4216295603543240908[14] = 1;
   out_4216295603543240908[15] = 0;
   out_4216295603543240908[16] = 0;
   out_4216295603543240908[17] = 0;
}
void h_30(double *state, double *unused, double *out_6673351104037270817) {
   out_6673351104037270817[0] = state[4];
}
void H_30(double *state, double *unused, double *out_9116354278712919174) {
   out_9116354278712919174[0] = 0;
   out_9116354278712919174[1] = 0;
   out_9116354278712919174[2] = 0;
   out_9116354278712919174[3] = 0;
   out_9116354278712919174[4] = 1;
   out_9116354278712919174[5] = 0;
   out_9116354278712919174[6] = 0;
   out_9116354278712919174[7] = 0;
   out_9116354278712919174[8] = 0;
}
void h_26(double *state, double *unused, double *out_6288308938331448344) {
   out_6288308938331448344[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3070553517615327591) {
   out_3070553517615327591[0] = 0;
   out_3070553517615327591[1] = 0;
   out_3070553517615327591[2] = 0;
   out_3070553517615327591[3] = 0;
   out_3070553517615327591[4] = 0;
   out_3070553517615327591[5] = 0;
   out_3070553517615327591[6] = 0;
   out_3070553517615327591[7] = 1;
   out_3070553517615327591[8] = 0;
}
void h_27(double *state, double *unused, double *out_6170301544398071600) {
   out_6170301544398071600[0] = state[3];
}
void H_27(double *state, double *unused, double *out_6892760207528975957) {
   out_6892760207528975957[0] = 0;
   out_6892760207528975957[1] = 0;
   out_6892760207528975957[2] = 0;
   out_6892760207528975957[3] = 1;
   out_6892760207528975957[4] = 0;
   out_6892760207528975957[5] = 0;
   out_6892760207528975957[6] = 0;
   out_6892760207528975957[7] = 0;
   out_6892760207528975957[8] = 0;
}
void h_29(double *state, double *unused, double *out_6628502906556633262) {
   out_6628502906556633262[0] = state[1];
}
void H_29(double *state, double *unused, double *out_8606122934398526990) {
   out_8606122934398526990[0] = 0;
   out_8606122934398526990[1] = 1;
   out_8606122934398526990[2] = 0;
   out_8606122934398526990[3] = 0;
   out_8606122934398526990[4] = 0;
   out_8606122934398526990[5] = 0;
   out_8606122934398526990[6] = 0;
   out_8606122934398526990[7] = 0;
   out_8606122934398526990[8] = 0;
}
void h_28(double *state, double *unused, double *out_2850241868042404842) {
   out_2850241868042404842[0] = state[0];
}
void H_28(double *state, double *unused, double *out_4758222122241494052) {
   out_4758222122241494052[0] = 1;
   out_4758222122241494052[1] = 0;
   out_4758222122241494052[2] = 0;
   out_4758222122241494052[3] = 0;
   out_4758222122241494052[4] = 0;
   out_4758222122241494052[5] = 0;
   out_4758222122241494052[6] = 0;
   out_4758222122241494052[7] = 0;
   out_4758222122241494052[8] = 0;
}
void h_31(double *state, double *unused, double *out_3408709555238244804) {
   out_3408709555238244804[0] = state[8];
}
void H_31(double *state, double *unused, double *out_6842702798366344243) {
   out_6842702798366344243[0] = 0;
   out_6842702798366344243[1] = 0;
   out_6842702798366344243[2] = 0;
   out_6842702798366344243[3] = 0;
   out_6842702798366344243[4] = 0;
   out_6842702798366344243[5] = 0;
   out_6842702798366344243[6] = 0;
   out_6842702798366344243[7] = 0;
   out_6842702798366344243[8] = 1;
}
#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;

  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);

  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H);

  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();


    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;

  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);

  // update cov
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}




}
extern "C" {

void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
}
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<2, 3, 0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
}
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
}
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
}
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
}
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
}
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
}
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
  update<1, 3, 0>(in_x, in_P, h_31, H_31, NULL, in_z, in_R, in_ea, MAHA_THRESH_31);
}
void car_err_fun(double *nom_x, double *delta_x, double *out_5872706819683017341) {
  err_fun(nom_x, delta_x, out_5872706819683017341);
}
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4273870982427227725) {
  inv_err_fun(nom_x, true_x, out_4273870982427227725);
}
void car_H_mod_fun(double *state, double *out_1982852844921793486) {
  H_mod_fun(state, out_1982852844921793486);
}
void car_f_fun(double *state, double dt, double *out_8224871313240091039) {
  f_fun(state,  dt, out_8224871313240091039);
}
void car_F_fun(double *state, double dt, double *out_5540635882771833240) {
  F_fun(state,  dt, out_5540635882771833240);
}
void car_h_25(double *state, double *unused, double *out_3133515492953738915) {
  h_25(state, unused, out_3133515492953738915);
}
void car_H_25(double *state, double *unused, double *out_6812056836489383815) {
  H_25(state, unused, out_6812056836489383815);
}
void car_h_24(double *state, double *unused, double *out_9067730718024463655) {
  h_24(state, unused, out_9067730718024463655);
}
void car_H_24(double *state, double *unused, double *out_4216295603543240908) {
  H_24(state, unused, out_4216295603543240908);
}
void car_h_30(double *state, double *unused, double *out_6673351104037270817) {
  h_30(state, unused, out_6673351104037270817);
}
void car_H_30(double *state, double *unused, double *out_9116354278712919174) {
  H_30(state, unused, out_9116354278712919174);
}
void car_h_26(double *state, double *unused, double *out_6288308938331448344) {
  h_26(state, unused, out_6288308938331448344);
}
void car_H_26(double *state, double *unused, double *out_3070553517615327591) {
  H_26(state, unused, out_3070553517615327591);
}
void car_h_27(double *state, double *unused, double *out_6170301544398071600) {
  h_27(state, unused, out_6170301544398071600);
}
void car_H_27(double *state, double *unused, double *out_6892760207528975957) {
  H_27(state, unused, out_6892760207528975957);
}
void car_h_29(double *state, double *unused, double *out_6628502906556633262) {
  h_29(state, unused, out_6628502906556633262);
}
void car_H_29(double *state, double *unused, double *out_8606122934398526990) {
  H_29(state, unused, out_8606122934398526990);
}
void car_h_28(double *state, double *unused, double *out_2850241868042404842) {
  h_28(state, unused, out_2850241868042404842);
}
void car_H_28(double *state, double *unused, double *out_4758222122241494052) {
  H_28(state, unused, out_4758222122241494052);
}
void car_h_31(double *state, double *unused, double *out_3408709555238244804) {
  h_31(state, unused, out_3408709555238244804);
}
void car_H_31(double *state, double *unused, double *out_6842702798366344243) {
  H_31(state, unused, out_6842702798366344243);
}
void car_predict(double *in_x, double *in_P, double *in_Q, double dt) {
  predict(in_x, in_P, in_Q, dt);
}
void car_set_mass(double x) {
  set_mass(x);
}
void car_set_rotational_inertia(double x) {
  set_rotational_inertia(x);
}
void car_set_center_to_front(double x) {
  set_center_to_front(x);
}
void car_set_center_to_rear(double x) {
  set_center_to_rear(x);
}
void car_set_stiffness_front(double x) {
  set_stiffness_front(x);
}
void car_set_stiffness_rear(double x) {
  set_stiffness_rear(x);
}
}

const EKF car = {
  .name = "car",
  .kinds = { 25, 24, 30, 26, 27, 29, 28, 31 },
  .feature_kinds = {  },
  .f_fun = car_f_fun,
  .F_fun = car_F_fun,
  .err_fun = car_err_fun,
  .inv_err_fun = car_inv_err_fun,
  .H_mod_fun = car_H_mod_fun,
  .predict = car_predict,
  .hs = {
    { 25, car_h_25 },
    { 24, car_h_24 },
    { 30, car_h_30 },
    { 26, car_h_26 },
    { 27, car_h_27 },
    { 29, car_h_29 },
    { 28, car_h_28 },
    { 31, car_h_31 },
  },
  .Hs = {
    { 25, car_H_25 },
    { 24, car_H_24 },
    { 30, car_H_30 },
    { 26, car_H_26 },
    { 27, car_H_27 },
    { 29, car_H_29 },
    { 28, car_H_28 },
    { 31, car_H_31 },
  },
  .updates = {
    { 25, car_update_25 },
    { 24, car_update_24 },
    { 30, car_update_30 },
    { 26, car_update_26 },
    { 27, car_update_27 },
    { 29, car_update_29 },
    { 28, car_update_28 },
    { 31, car_update_31 },
  },
  .Hes = {
  },
  .sets = {
    { "mass", car_set_mass },
    { "rotational_inertia", car_set_rotational_inertia },
    { "center_to_front", car_set_center_to_front },
    { "center_to_rear", car_set_center_to_rear },
    { "stiffness_front", car_set_stiffness_front },
    { "stiffness_rear", car_set_stiffness_rear },
  },
  .extra_routines = {
  },
};

ekf_init(car);
