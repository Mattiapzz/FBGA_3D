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

#ifndef NODESTRUCT3D_HXX
#define NODESTRUCT3D_HXX

#include "types.hxx"
#include <cmath>


namespace GG
{

  template<typename T>
  constexpr T eval_v(const T& S, const T& A, const T& V0) {
      return std::sqrt(static_cast<T>(2) * A * S + std::pow(V0, 2));
  }

  template<typename T>
  constexpr T eval_v2(const T& S, const T& A, const T& V0) {
      return static_cast<T>(2) * A * S + std::pow(V0, 2);
  }

  struct NodeStruct3D
  {
    // STATIC MEMBERS (GIVEN)
    //// Length 
    real s{0};
    //// geometry (given)
    ///// euler angles
    real mu{0.0};
    real phi{0.0};
    real theta{0.0};
    ///// euler derivatives
    real mu_prime{0.0};
    real phi_prime{0.0};
    real theta_prime{0.0};
    ///// elurer secodn derivatives
    real mu_double_prime{0.0};
    real phi_double_prime{0.0};
    real theta_double_prime{0.0};
    ////// offset
    real n{0.0};
    real chi{0.0};
    ////// additional offsets
    real chi_prime{0.0};
    // STATIC MEMBERS (COMPUTED)
    //// Gravity corrections
    real g_x{0.0};
    real g_y{0.0};
    real g_z{0.0};
    //// Geometric Omegas 
    real Omega_x{0.0};
    real Omega_y{0.0};
    real Omega_z{0.0};
    //// Geometric Omegas prime
    real Omega_x_prime{0.0};
    real Omega_y_prime{0.0};
    real Omega_z_prime{0.0};
    //
    real V_max{130.0};
    // Adherence scaling
    real alpha{1.0};
  };

  struct CellStruct3D
  {
    real s_0{0};
    real s_1{0};
    real L{0};
    size_type ID0{0};
    size_type ID1{0};
    SegmentType m_type{UNKNOWN};
    real V_max{130.0};
    real V_dot{QUIET_NAN};
    real V_0{0.0};
    real V_1{0.0};
  };

}


#endif // NODESTRUCT3D_HXX