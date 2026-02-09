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

#pragma once

#include <vector>

#include <cstddef>
#include <limits>
#ifndef TYPES_HXX
#define TYPES_HXX

#define STD_TOL 1.0e-10
#define STD_MAX_ITER 200
#define STD_VERBOSE "zero"

namespace GG
{

  /*\
   |   _____                     _       __
   |  |_   _|   _ _ __   ___  __| | ___ / _|___
   |    | || | | | '_ \ / _ \/ _` |/ _ \ |_/ __|
   |    | || |_| | |_) |  __/ (_| |  __/  _\__ \
   |    |_| \__, | .__/ \___|\__,_|\___|_| |___/
   |        |___/|_|
  \*/

  // Enum to understand if a segment is of type forward (1) or backward (2) or transitions (3) or unknown (0)
  enum SegmentType
  {
    FORWARD = 1,
    BACKWARD = 2,
    TRANSITION = 3,
    FORWARD_NAN = 4,
    BACKWARD_NAN = 5,
    TRANSITION_NAN = 6,
    YELLOWFLAG = 7,
    YELLOWFLAG_NAN = 8,
    UNKNOWN = 0
  };

  using real = double;             //!< Real number type
  using floating = float;          //!< Real number type
  using integer = int;             //!< Integer number type
  using size_type = std::size_t;   //!< Size type

  /*\
   |    ____                _              _
   |   / ___|___  _ __  ___| |_ __ _ _ __ | |_ ___
   |  | |   / _ \| '_ \/ __| __/ _` | '_ \| __/ __|
   |  | |__| (_) | | | \__ \ || (_| | | | | |_\__ \
   |   \____\___/|_| |_|___/\__\__,_|_| |_|\__|___/
   |
  \*/

  static constexpr real GRAVITY = 9.81;                                   //!< Gravity static constant value

  static constexpr real EPSILON_MACHINE = std::numeric_limits<real>::epsilon(); //!< Machine epsilon epsilon static constant value
  static constexpr real EPSILON_HIGH = 1.0E-16;                                 //!< High precision epsilon static constant value
  static constexpr real EPSILON_MEDIUM = 1.0E-10;                               //!< Medium precision epsilon static constant value
  static constexpr real EPSILON_LOW = 1.0E-07;                                  //!< Low precision epsilon static constant value
  static constexpr real EPSILON = EPSILON_MEDIUM;                               //!< Standard precision epsilon static constant value
  static real const INFTY = std::numeric_limits<real>::infinity();          //!< Infinity static constant value
  static real const QUIET_NAN = std::numeric_limits<real>::quiet_NaN();     //!< Not-a-Number static constant value
  static constexpr real PI = 3.141592653589793238462643383279500;         //!< Pi static constant value
  static constexpr real PIDIV180 = 0.017453292519943295769236907684886;   //!< Pi/180 static constant value
  static constexpr real DEG2RAD = PIDIV180;                               //!< Degrees to Gradians static constant value
  static constexpr real RAD2DEG = 1.0 / DEG2RAD;                          //!< Gradians to Degrees static constant value

  struct scaling_gggv_factors
  {
    real ax_max_scale = 1.0; // Scale factor for a_x_max
    real ax_min_scale = 1.0; // Scale factor for a_x_min
    real ay_scale = 1.0;     // Scale factor for a_y
    real gg_exponent_ax_pos = 1.3; // Exponent for a_x_max positive
    real gg_exponent_ax_neg = 1.3; // Exponent for a_x
  };

  struct spline_data_collection
  {
    std::vector<real> v_data{0.0,90.0};
    std::vector<real> az_data{5.0,15.0};
    std::vector<real> ay_max_data{10.0,10.0,10.0,10.0};
    std::vector<real> ax_min_data{-10.0,-10.0,-10.0,-10.0};
    std::vector<real> ax_max_data{10.0,10.0,10.0,10.0};
    scaling_gggv_factors scaling_factors;
    std::vector<GG::real> v_data_eng{0.0,90.0}; // Velocity data for engine max GGGV
    std::vector<GG::real> ax_eng_data{13.0, 0.0}; // Acceleration data for engine max GGGV
  };

  struct engine_max_ggv
  {
    std::vector<GG::real> v_data_eng{0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0}; // Velocity data for engine max GGGV
    std::vector<GG::real> ax_eng_data{5.0, 10.0, 13.0, 10.0, 7.0, 4.6, 2.5, 0.3, 0.0, 0.0}; // Acceleration data for engine max GGGV
  };

  // [ 5.0, 10.0, 13.0, 10.0, 7.0, 4.6, 2.5, 0.3, 0.0, 0.0]
  // [10.0, 10.0,  9.5,  9.0, 7.4, 4.9, 2.9, 0.6, 0.1, 0.0]


  // stuct to include vector about the road angles and their geometric derivatives
  struct road_angles_and_derivatives_container {
    std::vector<real> mu; // Pitch angle
    std::vector<real> phi;  // Roll angle
    std::vector<real> theta;   // Yaw angle
    std::vector<real> mu_prime; // Geometric derivative of the pitch angle
    std::vector<real> phi_prime; // Geometric derivative of the roll angle
    std::vector<real> theta_prime; // Geometric derivative of the yaw angle
    std::vector<real> abscissa; // Geometric second derivative of the pitch angle
  };

  struct trajectory_offset_container {
    std::vector<real> n; 
    std::vector<real> chi;
  };

  struct adherence_container {
    std::vector<real> alpha; // Adherence coefficient
  };


  struct trajectory_offset_and_angles_container {
    trajectory_offset_container offset;
    road_angles_and_derivatives_container reference;
    adherence_container adherence; // Adherence coefficients
  };

  struct solution_container {
    std::vector<real> V0; // Initial velocity
    std::vector<real> V1; // Final velocity
    std::vector<real> V_dot; // Velocity derivative
  };

  struct output_plot_ggv_shell {
    std::vector<floating> a_tilde_x_max;
    std::vector<floating> a_tilde_x_min;
    std::vector<floating> a_tilde_y;
    std::vector<floating> v;
  };

  struct input_plot_ggv_shell {
    std::vector<real> v = {0.0,10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0};
    real az = GRAVITY;
    real alpha = 1.0;
  };

  struct output_plot_ggv_traj {
    std::vector<floating> a_tilde_x;
    std::vector<floating> a_tilde_y;
    std::vector<floating> v;
  };

  constexpr real VMAX = 130.0; // Maximum velocity
  struct yellow_flag_data {
    real v_des_max = VMAX;
    real a_des_min = -100.0;
    bool is_yellow = false; // Flag to indicate if the yellow flag is active
  };

  struct ablation_flags{
    // TODO activate flags to compute the different quantities
    bool compute_chi_dot = true; // Flag to compute chi_dot
    bool compute_s_ddot = true; // Flag to compute s_ddot
    bool compute_Omega_xyz_prime = true; // Flag to compute Omega_xyz_prime
    bool compute_Euler_prime = true; // Flag to compute Euler angles derivatives
    bool compute_Euler_double_prime = true; // Flag to compute Euler angles second derivatives
  };

  struct constraint_violation{
    real maximum = -INFTY; // Maximum violation
    real minimum = +INFTY; // Minimum violation
    real average = QUIET_NAN; // Average violation
    integer id_max = -1; // ID of the maximum violation
    integer id_min = -1; // ID of the minimum violation
    real H_inf = INFTY; // H_inf norm of the violation
    real cumulative = 0.0; // Cumulative violation
    std::vector<real> values; // Vector of violations
    std::vector<real> values_violation;
    integer num_violations = 0; // Number of violations
  };

  struct constraint_violation_for_segments{
    std::vector<std::vector<real>> values;
    std::vector<std::vector<real>> values_violation;
    real maximum = -INFTY; // Maximum violation
    real minimum = INFTY; // Minimum violation
    real average = QUIET_NAN; // Average violation
    integer id_max = -1; // ID of the maximum violation
    integer id_min = -1; // ID of the minimum violation
    real H_inf = INFTY; // H_inf norm of the violation
    real cumulative = 0.0; // Cumulative violation
    std::vector<real> cumulative_by_segment; // Cumulative violation by segment
    std::vector<real> average_by_segment; // Average violation by segment
    std::vector<real> H_inf_by_segment; // H_inf norm of the violation by segment
    integer num_violations = 0; // Number of violations
  };

  struct Context {
    // Basic kinematic values
    real cos_chi{1}, sin_chi{0};
    real den_common{1}, inv_den_common{1}, inv_den_common2{1}, inv_den_common3{1};
    real s_dot{0}, w{0}, chi_dot{0}, n_dot{0};
    
    // Second derivatives
    real s_ddotA{0}, s_ddotB{0}, s_ddot{0}, w_dot{0};
    
    // Angular velocities
    real omega_hat_x{0}, omega_hat_y{0}, omega_hat_z{0};
    
    // Accelerations
    real a_tilde_x{0}, a_tilde_y{0}, a_tilde_z{0};
    
    // G-G diagram parameters
    real a_tilde_y_lim{0}, a_tilde_y_clip{0};
    real a_tilde_x_max_gg{0}, a_tilde_x_min_gg{0}, a_tilde_x_eng{0};
    real rho_max{0}, rho_min{0};
    
    // Final computed bounds
    real a_tilde_x_max{0}, a_tilde_x_min{0};
    
    // Edge case flag
    bool at_lateral_limit{0};
};

struct solver_params
{
  real tolerance        = STD_TOL;
  int max_iter          = STD_MAX_ITER;
  std::string verbosity = STD_VERBOSE;
};

} // namespace FBGA3D

#endif // TYPES_HXX