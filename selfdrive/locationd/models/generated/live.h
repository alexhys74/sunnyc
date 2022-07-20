#pragma once
#include "rednose/helpers/common_ekf.h"
extern "C" {
void live_update_4(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_9(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_10(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_12(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_31(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_32(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_13(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_14(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_update_33(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea);
void live_H(double *in_vec, double *out_4790583136941706646);
void live_err_fun(double *nom_x, double *delta_x, double *out_2909834052534932465);
void live_inv_err_fun(double *nom_x, double *true_x, double *out_8555278309126036397);
void live_H_mod_fun(double *state, double *out_3435397669119375336);
void live_f_fun(double *state, double dt, double *out_5039992466047587188);
void live_F_fun(double *state, double dt, double *out_6569067551416855978);
void live_h_4(double *state, double *unused, double *out_7957359290500711161);
void live_H_4(double *state, double *unused, double *out_3736089306800884482);
void live_h_9(double *state, double *unused, double *out_4030518270746111525);
void live_H_9(double *state, double *unused, double *out_3551129628463562988);
void live_h_10(double *state, double *unused, double *out_4457721155288270543);
void live_H_10(double *state, double *unused, double *out_5132166839472521157);
void live_h_12(double *state, double *unused, double *out_4228401283344149185);
void live_H_12(double *state, double *unused, double *out_8329396389865934138);
void live_h_31(double *state, double *unused, double *out_5208820582285086291);
void live_H_31(double *state, double *unused, double *out_369427249428277106);
void live_h_32(double *state, double *unused, double *out_4222745714023749077);
void live_H_32(double *state, double *unused, double *out_1571411087904731695);
void live_h_13(double *state, double *unused, double *out_4914072219343370484);
void live_H_13(double *state, double *unused, double *out_1121568688562748714);
void live_h_14(double *state, double *unused, double *out_4030518270746111525);
void live_H_14(double *state, double *unused, double *out_3551129628463562988);
void live_h_33(double *state, double *unused, double *out_8594229317744820794);
void live_H_33(double *state, double *unused, double *out_8619585029864114293);
void live_predict(double *in_x, double *in_P, double *in_Q, double dt);
}