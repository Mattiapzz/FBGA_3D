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

#ifndef FBGA_MOTO_HXX
#define FBGA_MOTO_HXX

#include "brentdekker.hxx"
#include "NodeStruct3D.hxx"
#include "types.hxx"
#include "FBGA3D/gggv_MOTO.hxx"
#include <vector>

#define DEBUG_METHODS 1

#define NUMPTSPLOTSHELL 100

namespace GG
{



#define CHECK_CONSTRAINTS_LEVEL 0 // 0: no final check, 1: final check on nodes, 2: final check between nodes
#define NUMPTS_CHECK_CONSTRAINTS 100 // Number of points to check constraints


class FBGA_MOTO
{
private:
  /* data */
  GG::gggv_MOTO gggv_moto; // GGv object
  brentdekker BD;                // Solver
  // minimizer MGS;              // Minimizer
  solver_params solver_p;        // Solver parameters
  std::vector<NodeStruct3D> Nodes; // Vector of NodeStruct3D
  std::vector<CellStruct3D> Cells; // Vector of CellStruct3D
  real v_I{0.0};                 // Initial velocity
  real v_max{VMAX};              // Maximum velocity
  std::vector<size_type> dump_seg_id;  // Vector of segments with problems for debug
  integer m_past_index = 0;
  real m_lat_tol = 1e-6; //1e-6;
  ablation_flags m_ablation_flags; // Ablation flags
  real m_lateral_shrink_factor = 1.0;
  real m_lat_tol_vmax = 1e-4; //1e-6;
  //
  integer yellow_index = 0;
  yellow_flag_data m_yellow_flag_data; // Yellow flag data

  #if CHECK_CONSTRAINTS_LEVEL > 0
  constraint_violation m_constraint_violation; // Constraint violation
  #if CHECK_CONSTRAINTS_LEVEL > 1
  constraint_violation_for_segments m_constraint_violation_segments; // Constraint violation for segments
  #endif
  #endif


public:
  // constructors
  FBGA_MOTO();

  real compute(trajectory_offset_and_angles_container const &TOA, real V0 = -1.0, yellow_flag_data const & yellow_flag = yellow_flag_data());

private:
  // Setup method
  void create_nodes_cells(trajectory_offset_and_angles_container const &TOA);
  // compute Vmax vector
  void compute_Vmax();
  // FW
  void FW();
  // BW
  void BW();
  // BY
  void BY();

  #if CHECK_CONSTRAINTS_LEVEL > 0

  void final_constraints_check();
  
  #if CHECK_CONSTRAINTS_LEVEL > 1
  void final_constraints_check_segments();
  #endif // CHECK_CONSTRAINTS_LEVEL > 1
  #endif // CHECK_CONSTRAINTS_LEVEL > 0


  //
  [[nodiscard]] real compute_total_time() const;

  [[nodiscard]] real compute_time(real s) const;


  // static method to evaluate Omega_xyz 
  static void eval_Omega_xyz(GG::NodeStruct3D & node);
  // static method to evaluate Omega_prime_xyz
  static void eval_Omega_prime_xyz(GG::NodeStruct3D & node);
  // static method to evaluate g_xyz
  static void eval_g_xyz(GG::NodeStruct3D & node);

  static real eval_a_tilde_x(GG::NodeStruct3D const & node, real V, real V_dot);
  static real eval_a_hat_x(GG::NodeStruct3D const & node, real V, real V_dot);
  static real eval_V_dot_Vatildex(GG::NodeStruct3D const & node, real V, real atildex);
  static real eval_V_dot_Vahatx(GG::NodeStruct3D const & node, real V, real ahatx);
  
  [[nodiscard]] real eval_a_tilde_y(GG::NodeStruct3D const & node, real V);
  [[nodiscard]] real eval_a_hat_y(GG::NodeStruct3D const & node, real V);
  
  [[nodiscard]] real eval_a_tilde_z(GG::NodeStruct3D const & node, real V, real V_dot);
  [[nodiscard]] real eval_a_hat_z(GG::NodeStruct3D const & node, real V, real V_dot);
  
  
  
  [[nodiscard]] real eval_a_tilde_x_max(GG::NodeStruct3D const & node, real V, real V_dot);
  [[nodiscard]] real eval_a_tilde_x_min(GG::NodeStruct3D const & node, real V, real V_dot);
  [[nodiscard]] real eval_a_tilde_y_lim(GG::NodeStruct3D const & node, real V, real V_dot);

  [[nodiscard]] real eval_a_tilde_x_min_YF(GG::NodeStruct3D const & node, real V, real V_dot);


  [[nodiscard]]real eval_s_ddot(GG::NodeStruct3D const & node, real V, real V_dot);


  [[nodiscard]] real signed_distance(GG::NodeStruct3D const & node, real V, real V_dot);
  [[nodiscard]] real signed_distance_YF(GG::NodeStruct3D const & node, real V, real V_dot);

  [[nodiscard]] real signed_distanceNOOPTI(GG::NodeStruct3D const & node, real V, real V_dot);

  [[nodiscard]] static real eval_V_next(GG::NodeStruct3D const & node, GG::CellStruct3D const & cell, real VDOT);
  [[nodiscard]] static real eval_V_prev(GG::NodeStruct3D const & node, GG::CellStruct3D const & cell, real VDOT);
  [[nodiscard]] static real eval_V(GG::NodeStruct3D const & node, GG::CellStruct3D const & cell, real s);

  [[nodiscard]] real max_lateral_performance_func(GG::NodeStruct3D const & node, real V);

  [[nodiscard]]real max_longitudinal_performance_func(GG::NodeStruct3D const & node, real V, real V_dot);
  [[nodiscard]]real max_longitudinal_performance_funcNOOPTI(GG::NodeStruct3D const & node, real V, real V_dot);

  [[nodiscard]]real get_max_longitudinal_performance(GG::NodeStruct3D const & node, real V);

  [[nodiscard]]real min_longitudinal_performance_func(GG::NodeStruct3D const & node, real V, real V_dot);
  [[nodiscard]]real min_longitudinal_performance_funcNOOPTI(GG::NodeStruct3D const & node, real V, real V_dot);

  [[nodiscard]]real get_min_longitudinal_performance(GG::NodeStruct3D const & node, real V);


  [[nodiscard]]real min_longitudinal_performance_func_YF(GG::NodeStruct3D const & node, real V, real V_dot);

  [[nodiscard]]real get_min_longitudinal_performance_YF(GG::NodeStruct3D const & node, real V);

  inline void compute_context(const GG::NodeStruct3D& node, real V, real V_dot, Context & ctx) const;


  [[nodiscard]] real eval_contraint_satisfaction(GG::NodeStruct3D const & node, real V, real V_dot);

public:

  void eval_contraint_satisfaction();

  [[nodiscard]] integer get_cell_idx(const real s) const;
  [[nodiscard]] real eval_Vmax(real s) const;
  [[nodiscard]] real eval_V(real s) const;
  [[nodiscard]] real eval_V_dot(real s) const;
  [[nodiscard]] real eval_A_hat_x(real s) const;
  [[nodiscard]] real eval_A_hat_y(real s) ;
  [[nodiscard]] real eval_A_hat_z(real s) ;  

  [[nodiscard]] real eval_A_tilde_x(real s) const;
  [[nodiscard]] real eval_A_tilde_y(real s) ;
  [[nodiscard]] real eval_A_tilde_z(real s) ;

  [[nodiscard]] real eval_g_x(real s) const;
  [[nodiscard]] real eval_g_y(real s) const;
  [[nodiscard]] real eval_g_z(real s) const;


  [[nodiscard]] real eval_Omega_x(real s) const;
  [[nodiscard]] real eval_Omega_y(real s) const;
  [[nodiscard]] real eval_Omega_z(real s) const;

  [[nodiscard]] SegmentType eval_segment_type(real s) const;
  
  [[nodiscard]] real eval_alpha(real s) const;


  [[nodiscard]] output_plot_ggv_shell eval_shell_plot(GG::NodeStruct3D const & fictitious_node, real Va, real Vb, real V_dot) ;

  [[nodiscard]] output_plot_ggv_shell eval_shell_plotpy(GG::input_plot_ggv_shell const & input_plot) ;


  // static method to evaluate Omega_xyz 
  void eval_Omega_xyz_plot(GG::NodeStruct3D & fictitious_node) const
  {
    this->eval_Omega_xyz(fictitious_node);
  }
  // // static method to evaluate g_xyz
  void eval_g_xyz_plot(GG::NodeStruct3D & fictitious_node) const
  {
    this->eval_g_xyz(fictitious_node);
  }



};



} // namespace GG

#endif // FBGA_MOTO_HXX



