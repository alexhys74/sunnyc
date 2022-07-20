#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void car_update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void car_err_fun(double *nom_x, double *delta_x, double *out_5872706819683017341);
void car_inv_err_fun(double *nom_x, double *true_x, double *out_4273870982427227725);
void car_H_mod_fun(double *state, double *out_1982852844921793486);
void car_f_fun(double *state, double dt, double *out_8224871313240091039);
void car_F_fun(double *state, double dt, double *out_5540635882771833240);
void car_h_25(double *state, double *unused, double *out_3133515492953738915);
void car_H_25(double *state, double *unused, double *out_6812056836489383815);
void car_h_24(double *state, double *unused, double *out_9067730718024463655);
void car_H_24(double *state, double *unused, double *out_4216295603543240908);
void car_h_30(double *state, double *unused, double *out_6673351104037270817);
void car_H_30(double *state, double *unused, double *out_9116354278712919174);
void car_h_26(double *state, double *unused, double *out_6288308938331448344);
void car_H_26(double *state, double *unused, double *out_3070553517615327591);
void car_h_27(double *state, double *unused, double *out_6170301544398071600);
void car_H_27(double *state, double *unused, double *out_6892760207528975957);
void car_h_29(double *state, double *unused, double *out_6628502906556633262);
void car_H_29(double *state, double *unused, double *out_8606122934398526990);
void car_h_28(double *state, double *unused, double *out_2850241868042404842);
void car_H_28(double *state, double *unused, double *out_4758222122241494052);
void car_h_31(double *state, double *unused, double *out_3408709555238244804);
void car_H_31(double *state, double *unused, double *out_6842702798366344243);
void car_predict(double *in_x, double *in_P, double *in_Q, double dt);
void car_set_mass(double x);
void car_set_rotational_inertia(double x);
void car_set_center_to_front(double x);
void car_set_center_to_rear(double x);
void car_set_stiffness_front(double x);
void car_set_stiffness_rear(double x);
}