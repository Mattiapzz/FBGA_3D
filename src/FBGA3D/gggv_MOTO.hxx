/*
(***********************************************************************************)
(*                                                                                 *)
(* The FBGA3D project                                                              *)
(*                                                                                 *)
(* Copyright (c) 2025, Mattia Piazza                                               *)
(*                                                                                 *)
(* The FBGA3D project and its components are supplied under the terms              *)
(* of the open source BSD 2-Clause License. The contents of the FBGA3D             *)
(* project and its components may not be copied or disclosed except in             *)
(* accordance with the terms of the BSD 2-Clause License.                          *)
(*                                                                                 *)
(* URL: https://opensource.org/licenses/BSD-2-Clause                               *)
(*                                                                                 *)
(*    Mattia Piazza                                                                *)
(*    Department of Industrial Engineering                                         *)
(*    University of Trento                                                         *)
(*    e-mail: mattia.piazza@unitn.it                                               *)
(*                                                                                 *)
(***********************************************************************************)
*/

// PerformanceEnvelope.h
#ifndef GGGV_MOTO_HXX
#define GGGV_MOTO_HXX


#include "types.hxx"

namespace GG {

  // data structure to hold the aerodynamic data
  struct MotoData{
    GG::real b = 0.73; // wheelbase
    GG::real L_W = 1.5; // length of the wheel
    GG::real h = 0.69; // height of the center of mass
    GG::real mu_X = 1.30; // longitudinal friction coefficient
    GG::real mu_Y = 1.40; // lateral friction coefficient
    GG::real c_a_0 = 0.00; // aerodynamic coefficient
    GG::real c_a_1 = 0.00; // aerodynamic coefficient
    GG::real c_a_2 = 0.0*0.5*1.2*0.25/(250.0*9.81); // aerodynamic coefficient
    GG::real h_a = 0.51; // height of the aerodynamic center
    GG::real g = 9.81; // gravitational acceleration
    GG::real M = 250.0; // mass of the motorcycle
    GG::real P = 145.0 * 1000.0; // maximum power of the engine in Watts
  };

  class gggv_MOTO { 
    private:


    MotoData m_aero_data; // Aerodynamic data
    // real m_g = GRAVITY;

    public:
    gggv_MOTO();
    // Pure virtual method to be overridden by derived classes
    [[nodiscard]] real a_x_push( real ay_tilde, real V, real az_tilde) const;
    [[nodiscard]] real a_x_pull( real ay_tilde, real V, real az_tilde) const;

    [[nodiscard]] real a_x_eng(real V) const;
    [[nodiscard]] real a_y_lim(real V, real az_tilde) const;

    [[nodiscard]] real a_x_aero( real V) const;
    
    //
    // private:
    void setup_std();

    [[nodiscard]] real constraint_ax_max_wheeling(real ay_hat, real az_hat, real v) const;
    [[nodiscard]] real constraint_ax_max_adherence(real ay_hat, real az_hat, real v) const;
    [[nodiscard]] real constraint_ax_min_stoppie(real ay_hat, real az_hat, real v) const;
    [[nodiscard]] real constraint_ax_min_adherence(real ay_hat, real az_hat, real v) const;
    [[nodiscard]] real constraint_ay_lim_adherence(real az_hat, real v) const;

  };
  
}
#endif // GGGV_MOTO_HXX