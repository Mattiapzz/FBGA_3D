#include "FBGA3D/FBGA_MOTO.hxx"
#include "FBGA3D/gg_utils.hxx"
#include "FBGA3D/NodeStruct3D.hxx"
#include "FBGA3D/types.hxx"
#include <algorithm>
#include <cmath>
#include <cstddef>

#define DEBUG_FB 0

#include <iostream>

using namespace GG;

constexpr size_t DEFAULT_SIZE{100};

// --------------------------------------------------------------------------------------------

FBGA_MOTO::FBGA_MOTO()
{
  this->Nodes.reserve(DEFAULT_SIZE);
  this->Cells.reserve(DEFAULT_SIZE);
  this->dump_seg_id.reserve(DEFAULT_SIZE);
}

// --------------------------------------------------------------------------------------------

//  ██████╗███████╗██╗     ██╗                 ███╗   ██╗ ██████╗ ██████╗ ███████╗███████╗
// ██╔════╝██╔════╝██║     ██║                 ████╗  ██║██╔═══██╗██╔══██╗██╔════╝██╔════╝
// ██║     █████╗  ██║     ██║      █████╗     ██╔██╗ ██║██║   ██║██║  ██║█████╗  ███████╗
// ██║     ██╔══╝  ██║     ██║      ╚════╝     ██║╚██╗██║██║   ██║██║  ██║██╔══╝  ╚════██║
// ╚██████╗███████╗███████╗███████╗            ██║ ╚████║╚██████╔╝██████╔╝███████╗███████║
//  ╚═════╝╚══════╝╚══════╝╚══════╝            ╚═╝  ╚═══╝ ╚═════╝ ╚═════╝ ╚══════╝╚══════╝
                                                                                  
// --------------------------------------------------------------------------------------------

void 
FBGA_MOTO::create_nodes_cells(trajectory_offset_and_angles_container const &TOA)
{
  this->Nodes.clear();
  this->Cells.clear();  
  const size_type N = TOA.reference.abscissa.size();
  this->Nodes.resize(N);
  this->Cells.resize(N-1);

  std::vector<real> tmp_dchi         = computeFiniteDifference(TOA.reference.abscissa, TOA.offset.chi);
  std::vector<real> tmp_dmu_prime    = computeFiniteDifference(TOA.reference.abscissa, TOA.reference.mu_prime);
  std::vector<real> tmp_dphi_prime   = computeFiniteDifference(TOA.reference.abscissa, TOA.reference.phi_prime);
  std::vector<real> tmp_dtheta_prime = computeFiniteDifference(TOA.reference.abscissa, TOA.reference.theta_prime);
  
  //
  for(size_type i = 0; i < (N); i++)
  {
    auto & node = this->Nodes[i];
    //
    node.s           = TOA.reference.abscissa[i];
    //
    node.mu          = TOA.reference.mu[i];
    node.phi         = TOA.reference.phi[i];
    node.theta       = TOA.reference.theta[i];
    //
    if(this->m_ablation_flags.compute_Euler_prime) // do not compute in ablation analysis mode
    {
      node.mu_prime    = TOA.reference.mu_prime[i];
      node.phi_prime   = TOA.reference.phi_prime[i];
      node.theta_prime = TOA.reference.theta_prime[i];
    }
    // //
    if(this->m_ablation_flags.compute_Euler_double_prime) // do not compute in ablation analysis mode
    {
      node.mu_double_prime    = tmp_dmu_prime[i];
      node.phi_double_prime   = tmp_dphi_prime[i];
      node.theta_double_prime = tmp_dtheta_prime[i];
    }
    //
    node.n           = TOA.offset.n[i];
    node.chi         = TOA.offset.chi[i];
    //
    node.chi_prime   = tmp_dchi[i];
    //
    node.alpha       = TOA.adherence.alpha[i];
    //
    this->eval_Omega_xyz(node);
    this->eval_g_xyz(node);
    if(this->m_ablation_flags.compute_Omega_xyz_prime) // do not compute in ablation analysis mode
    {
      this->eval_Omega_prime_xyz(node);
    }
    if( i > 0)
    {
      auto & cell = this->Cells[i-1];
      //
      cell.s_0 = TOA.reference.abscissa[i-1];
      cell.s_1 = TOA.reference.abscissa[i];
      cell.L   = cell.s_1 - cell.s_0;
      //
      cell.ID0 = i-1;
      cell.ID1 = i;
      //
      cell.m_type = SegmentType::UNKNOWN; 
      //
    }
  }
}

// --------------------------------------------------------------------------------------------

void 
FBGA_MOTO::eval_Omega_xyz(GG::NodeStruct3D & node)
{
  // Cache trigonometric values
  const real cos_mu  = std::cos(node.mu);
  const real sin_mu  = std::sin(node.mu);
  const real cos_phi = std::cos(node.phi);
  const real sin_phi = std::sin(node.phi);
  // compute the angular velocity vector
  node.Omega_x = node.phi_prime  - sin_mu * node.theta_prime;
  node.Omega_y = cos_mu * node.mu_prime + cos_mu  * sin_phi * node.theta_prime;
  node.Omega_z = -sin_phi * node.mu_prime + cos_mu  * cos_phi * node.theta_prime;
}

// --------------------------------------------------------------------------------------------

void 
FBGA_MOTO::eval_Omega_prime_xyz(GG::NodeStruct3D & node)
{
  
  // Cache trigonometric values
  const real cos_mu  = std::cos(node.mu);
  const real sin_mu  = std::sin(node.mu);
  const real cos_phi = std::cos(node.phi);
  const real sin_phi = std::sin(node.phi);
  // compute the angular velocity vector
  node.Omega_x_prime = node.phi_double_prime  - cos_mu * node.mu_prime * node.theta_prime - sin_mu * node.theta_double_prime;
  const real t1 = cos_mu * node.theta_double_prime - node.mu_prime * (sin_mu * node.theta_prime + node.phi_prime);
  const real t2 = cos_mu * node.theta_prime * node.phi_prime + node.mu_double_prime;
  node.Omega_y_prime = t1 * sin_phi + cos_phi * t2;
  node.Omega_z_prime = t1 * cos_phi - sin_phi * t2;
}

// --------------------------------------------------------------------------------------------

void
FBGA_MOTO::eval_g_xyz(GG::NodeStruct3D & node)
{
  // compute the gravity vector
  const real cos_mu  = std::cos(node.mu);
  const real sin_mu  = std::sin(node.mu);
  const real cos_phi = std::cos(node.phi);
  const real sin_phi = std::sin(node.phi);
  const real cos_chi = std::cos(node.chi);
  const real sin_chi = std::sin(node.chi);
  node.g_x = GRAVITY * (-sin_mu * cos_chi + cos_mu * sin_phi * sin_chi);
  node.g_y = GRAVITY * (sin_mu * sin_chi + cos_mu * sin_phi * cos_chi);
  node.g_z = GRAVITY * cos_mu * cos_phi;
}

// --------------------------------------------------------------------------------------------

//  ██████╗ ██████╗ ██████╗ ███████╗
// ██╔════╝██╔═══██╗██╔══██╗██╔════╝
// ██║     ██║   ██║██████╔╝█████╗  
// ██║     ██║   ██║██╔══██╗██╔══╝  
// ╚██████╗╚██████╔╝██║  ██║███████╗
//  ╚═════╝ ╚═════╝ ╚═╝  ╚═╝╚══════╝
                                 
// --------------------------------------------------------------------------------------------

real FBGA_MOTO::compute(trajectory_offset_and_angles_container const &TOA, const real v0, yellow_flag_data const & yellow_flag )
{
  this->m_yellow_flag_data = yellow_flag;
  this->dump_seg_id.clear();
  this->v_I = v0;
  //
  if (this->v_I <= 0.0)
  {
    this->v_I = this->v_max;
  }
  //
  this->create_nodes_cells(TOA);
  //
  this->compute_Vmax();
  //

  this->yellow_index = 0;
  if (this->m_yellow_flag_data.is_yellow)
  {
    // find the braking point last and update the values of VMAX accordingly
    this->BY();
  }
  // this->v_max = v_ref;
  //
  this->FW();
  //
  this->BW();
  //
  return this->compute_total_time(); 
}

// --------------------------------------------------------------------------------------------

void FBGA_MOTO::compute_Vmax()
{
  // extimate a small value of curvature.
  constexpr real k_small     = 1e-6;
  constexpr real v_top_speed = 130;
  real vmax                  = QUIET_NAN;
  // chose bracketing interval
  // loop all nodes
  for(auto & node : this->Nodes)
  {
    if(std::abs(node.Omega_z)<k_small)
    {
      vmax = v_top_speed;
    }
    else
    {
      constexpr real v_a = 0.0;
      constexpr real v_b = v_top_speed;
      auto F2solve = [this,&node](const real v) -> real {
        return this->max_lateral_performance_func(node, v);
      };
      const bool ok = this->BD.solve(F2solve, v_a, v_b, vmax);
      //vmax              = std::isnan(vmax) ? v_top_speed : vmax;
      vmax              = ok ? vmax : v_top_speed;
    }
    node.V_max = vmax;
  }
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::max_lateral_performance_func(GG::NodeStruct3D const & node, real V)
{
  const real a_tilde_aereo = - this->gggv_moto.a_x_aero(V);
  const real V_dot = this->eval_V_dot_Vatildex(node, V, a_tilde_aereo);
  return std::abs(this->eval_a_tilde_y(node, V)) - (this->eval_a_tilde_y_lim(node,V, V_dot) + this->m_lat_tol - this->m_lat_tol_vmax )*(this->m_lateral_shrink_factor);
  // TODO: Optimize (look to previous commits)
}

// --------------------------------------------------------------------------------------------


real
FBGA_MOTO::get_max_longitudinal_performance(
  const GG::NodeStruct3D & node,
  const real V)
{

  // create a lambda function only of V_dot
  auto F2solve = [this, &node, V](const real V_dot) -> real {
    return this->max_longitudinal_performance_func(node, V, V_dot);
  };

  const real a_tilde_aereo = - this->gggv_moto.a_x_aero(V);
  const real V_dot_0 = this->eval_V_dot_Vatildex(node, V, a_tilde_aereo);
  const real a_tilde_x_top = this->eval_a_tilde_x_max(node, V, V_dot_0);
  real V_dot_max = this->eval_V_dot_Vatildex(node, V, a_tilde_x_top);
  //
  if ( (std::abs(V_dot_max-V_dot_0) <= this->solver_p.tolerance) || ( std::abs(F2solve(V_dot_max)) <= this->solver_p.tolerance) )
  {
    return V_dot_max;
  }
  if(std::abs(F2solve(V_dot_0)) <= this->solver_p.tolerance)
  {
    // if the function is already close to zero, return the value
    return V_dot_0;
  }
  // TODO: better trick to avoid the while loop
  while (F2solve(V_dot_max) < -this->solver_p.tolerance)
  {
    constexpr real scale_gain = 1.5;
    V_dot_max = (V_dot_max - V_dot_0)*scale_gain + V_dot_0;
  }
  real V_dot_sol = V_dot_0;
  auto ok = this->BD.solve(
    F2solve, 
    V_dot_0, // lower bound 
    V_dot_max,  // upper bound
    V_dot_sol
  );
  if(!ok)
  {
    std::cout << "FBGA_MOTO::get_max_longitudinal_performance: Error in solving the maximum longitudinal performance function for V = " 
              << V << ", a_tilde_x_top = " << a_tilde_x_top << "\n";
    std::cout << "\tF2solve(V_dot_0)   = " << F2solve(V_dot_0) << "\n";
    std::cout << "\tF2solve(V_dot_max) = " << F2solve(V_dot_max) << "\n";
    std::cout << "\tV_dot_sol          = " << V_dot_sol << "\n";
    std::cout << "\tV_dot_0            = " << V_dot_0 << "\n";
    std::cout << "\tV_dot_max          = " << V_dot_max << "\n";
  }
  return ok ? V_dot_sol : V_dot_0; // return the solution if found, otherwise return the initial value
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::get_min_longitudinal_performance(
  const GG::NodeStruct3D & node,
  const real V)
{
  // create a lambda function only of V_dot
  auto F2solve = [this, &node, V](const real V_dot) -> real {
    return this->min_longitudinal_performance_func(node, V, V_dot);
  };
  const real V_dot_0 = this->eval_V_dot_Vatildex(node, V, 0.0);
  const real a_tilde_x_bottom = this->eval_a_tilde_x_min(node, V, 0.0);
  real V_dot_min = this->eval_V_dot_Vatildex(node, V, a_tilde_x_bottom);
  if (std::abs(V_dot_min-V_dot_0) <= this->solver_p.tolerance)
  {
    return V_dot_min;
  }
  if (std::abs(F2solve(V_dot_min)) <= this->solver_p.tolerance)
  {
    // if the function is already close to zero, return the value
    return V_dot_min;
  }
  if (std::abs(F2solve(V_dot_0)) <= this->solver_p.tolerance)
  {
    // if the function is already close to zero, return the value
    return V_dot_0;
  }
  // TODO: better trick to avoid the while loop
  while (F2solve(V_dot_min) > +this->solver_p.tolerance)
  {
    // std::cout << "FBGA_MOTO::get_min_longitudinal_performance: Warning: V_dot_min = " 
    //           << V_dot_min << " is not a valid solution for V = " << V 
    //           << ", trying to increase it.\n";
    V_dot_min = (V_dot_min - V_dot_0)*1.5 + V_dot_0;
  }
  real V_dot_sol = V_dot_0;
  auto ok = this->BD.solve(
    F2solve, 
    V_dot_min, // lower bound 
    V_dot_0,  // upper bound
    V_dot_sol
  );
  if(!ok)
  {
    std::cout << "FBGA_MOTO::get_min_longitudinal_performance: Error in solving the minimum longitudinal performance function for V = " 
              << V << ", a_tilde_x_bottom = " << a_tilde_x_bottom << "\n";
    std::cout << "\tF2solve(V_dot_0)   = " << F2solve(V_dot_0) << "\n";
    std::cout << "\tF2solve(V_dot_min) = " << F2solve(V_dot_min) << "\n";
    std::cout << "\tV_dot_sol          = " << V_dot_sol << "\n";
    std::cout << "\tV_dot_0            = " << V_dot_0 << "\n";
    std::cout << "\tV_dot_min          = " << V_dot_min << "\n";
  }
  return ok ? V_dot_sol : V_dot_0; // return the solution if found, otherwise return the initial value
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::min_longitudinal_performance_funcNOOPTI(
  const GG::NodeStruct3D & node,
  const real V,
  const real V_dot)
{
  return this->eval_a_tilde_x(node, V, V_dot) - this->eval_a_tilde_x_min(node, V, V_dot) - this->solver_p.tolerance;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::min_longitudinal_performance_func(
  const GG::NodeStruct3D & node,
  const real V,
  const real V_dot)
{
  Context ctx;
  this->compute_context(node, V, V_dot, ctx);
  return ctx.a_tilde_x - ctx.a_tilde_x_min - this->solver_p.tolerance;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::get_min_longitudinal_performance_YF(
  const GG::NodeStruct3D & node,
  const real V)
{
  // create a lambda function only of V_dot
  auto F2solve = [this, &node, V](const real V_dot) -> real {
    return this->min_longitudinal_performance_func_YF(node, V, V_dot);
  };
  const real V_dot_0 = this->eval_V_dot_Vatildex(node, V, 0.0);
  const real a_tilde_x_bottom = this->eval_a_tilde_x_min_YF(node, V, 0.0);
  real V_dot_min = this->eval_V_dot_Vatildex(node, V, a_tilde_x_bottom);
  if (std::abs(V_dot_min-V_dot_0) <= this->solver_p.tolerance)
  {
    return V_dot_min;
  }
  if (std::abs(F2solve(V_dot_min)) <= this->solver_p.tolerance)
  {
    return V_dot_min;
  }
  if (std::abs(F2solve(V_dot_0)) <= this->solver_p.tolerance)
  {
    return V_dot_0;
  }
  while (F2solve(V_dot_min) > +this->solver_p.tolerance)
  {
    V_dot_min = (V_dot_min - V_dot_0)*1.5 + V_dot_0;
  }
  real V_dot_sol = V_dot_0;
  auto ok = this->BD.solve(
    F2solve, 
    V_dot_min, // lower bound 
    V_dot_0,  // upper bound
    V_dot_sol
  );
  if(!ok)
  {
    std::cout << "FBGA_MOTO::get_min_longitudinal_performance_YF: Error in solving the minimum longitudinal performance function for V = " 
              << V << ", a_tilde_x_bottom = " << a_tilde_x_bottom << "\n";
    std::cout << "\tF2solve(V_dot_0)   = " << F2solve(V_dot_0) << "\n";
    std::cout << "\tF2solve(V_dot_min) = " << F2solve(V_dot_min) << "\n";
    std::cout << "\tV_dot_sol          = " << V_dot_sol << "\n";
    std::cout << "\tV_dot_0            = " << V_dot_0 << "\n";
    std::cout << "\tV_dot_min          = " << V_dot_min << "\n";
  }
  return ok ? V_dot_sol : V_dot_0; // return the solution if found, otherwise return the initial value
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::min_longitudinal_performance_func_YF(
  const GG::NodeStruct3D & node,
  const real V,
  const real V_dot)
{
  return this->eval_a_tilde_x(node, V, V_dot) - this->eval_a_tilde_x_min_YF(node, V, V_dot) - this->solver_p.tolerance;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::max_longitudinal_performance_funcNOOPTI(
  const GG::NodeStruct3D & node,
  const real V,
  const real V_dot)
{
  return this->eval_a_tilde_x(node, V, V_dot) - this->eval_a_tilde_x_max(node, V, V_dot) + this->solver_p.tolerance;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::max_longitudinal_performance_func(
  const GG::NodeStruct3D & node,
  const real V,
  const real V_dot)
{
  Context ctx;
  this->compute_context(node, V, V_dot, ctx);
  return ctx.a_tilde_x - ctx.a_tilde_x_max + this->solver_p.tolerance;
}


// ------------------------------------------------------------------------------------------------

void
FBGA_MOTO::BY()
{
  bool ok_by = false;
  real v_0 = std::min(this->v_I, this->Nodes[0].V_max);
  // set the max velocity of first point to the initial velocity
  this->Nodes[0].V_max = v_0;
  // for(auto & cell : this->Cells)
  for(integer ith=0; ith < static_cast<integer>(this->Cells.size()); ith++)
  {
    auto & cell = this->Cells[ith];
    //
    auto & node_0 = this->Nodes[cell.ID0];
    auto & node_1 = this->Nodes[cell.ID1];
    cell.V_0 = v_0;
    //
    if (v_0 <= this->m_yellow_flag_data.v_des_max)
    {
      this->yellow_index = ith;
      ok_by = true;
      break;
    }
    //
    const real V_dot_0   = this->eval_V_dot_Vatildex(node_0, v_0, 0.0);
    //
    real V_dot_min = this->get_min_longitudinal_performance_YF(node_0, v_0);
    // compute V_dot_v_ref to get the deceleration needed to reach v_ref
    const real a_hat_x_v_ref = ((this->m_yellow_flag_data.v_des_max * this->m_yellow_flag_data.v_des_max) - (v_0 * v_0))/(2.0 * cell.L);
    const real V_dot_v_ref = this->eval_V_dot_Vahatx(node_0, v_0, a_hat_x_v_ref);
    //
    V_dot_min = std::min(V_dot_min,V_dot_0);
    V_dot_min = std::max(V_dot_min,V_dot_v_ref);
    //
    const real distance_V_dot_min_1 = this->signed_distance_YF(node_1, this->eval_V_next(node_0, cell, V_dot_min), V_dot_min);
    //
    real VDOTSOL = 0.0;
    if ((distance_V_dot_min_1 <= 0) )
    {
      VDOTSOL = V_dot_min;
    }
    else 
    {
      this->BD.solve(
        [this,&node_0, &node_1, &cell](const real VDOT) -> real {
          return this->signed_distance_YF(node_1, this->eval_V_next(node_0, cell, VDOT), VDOT) + this->solver_p.tolerance; },
          V_dot_min,
          V_dot_0,
          VDOTSOL );
    }
    // check if the solution is valid
    if (!std::isnan(VDOTSOL))
    {
      // update cell
      cell.V_dot = VDOTSOL;
      cell.V_1   = this->eval_V_next(node_0, cell, VDOTSOL);
      // update next v0 with the final velocity of the segment of the maximum reachable velocity
      v_0 = cell.V_1;
      // set type
      cell.m_type = SegmentType::YELLOWFLAG;
      // set Vmax of the following segment to be the one of the current cell
      node_1.V_max = std::min(node_1.V_max, cell.V_1);
    }
    else
    {
      // set type
      cell.m_type = SegmentType::YELLOWFLAG_NAN;
      this->dump_seg_id.push_back(cell.ID0);
      //set the next vmax
      node_1.V_max = std::min(node_1.V_max, v_0);
      // set v_0 to be the maximum reachable velocity
      v_0 = node_1.V_max;
      cell.V_1 = v_0;
      //
      // std::cout << "FBGA_MOTO::BY: Found a NaN. meaning that current YF requested acceleration is not enough!\n";
      this->yellow_index = ith;
      ok_by = false;
      break;
    }
  }
  // Pass all V max later than yellow_index to lower the maximum speed to ref speed
  for(integer ith=this->yellow_index; ith < static_cast<integer>(this->Cells.size()); ith++)
  {
    auto & cell = this->Cells[ith];
    //
    auto & node_1 = this->Nodes[cell.ID1];
    //
    // here we edit the max speed for all next values if ok_by is true. otherwise we should modify only those after the first time we encounter a value less tha v_ref
    if (node_1.V_max <= this->m_yellow_flag_data.v_des_max) {
    // Set flag to start limiting speeds from this point onwards
      ok_by = true;
    }
    if(ok_by)
    {
      node_1.V_max = std::min(node_1.V_max, this->m_yellow_flag_data.v_des_max);
    } 
  }
}

// --------------------------------------------------------------------------------------------

void
FBGA_MOTO::FW()
{
  real v_I = this->Cells[this->yellow_index].V_0;
  if (this->yellow_index==0)
  {
    v_I = this->v_I;
  }
  real v_0 = std::min(v_I, this->Nodes[this->yellow_index].V_max);
  // set the max velocity of first point to the initial velocity
  this->Nodes[this->yellow_index].V_max = v_0;
  // for(auto & cell : this->Cells)
  for(integer ith=this->yellow_index; ith < static_cast<integer>(this->Cells.size()); ith++)
  {
    auto & cell = this->Cells[ith];
    //
    auto & node_0 = this->Nodes[cell.ID0];
    auto & node_1 = this->Nodes[cell.ID1];
    cell.V_0 = v_0;
    // 
    const real a_hat_x_Vmax = ((node_1.V_max * node_1.V_max) - (v_0 * v_0))/(2.0 * cell.L);
    const real V_dot_x_Vmax = this->eval_V_dot_Vahatx(node_0, v_0, a_hat_x_Vmax);
    const real V_dot_0   = this->eval_V_dot_Vatildex(node_0, v_0, 0.0);
    //
    real V_dot_max = this->get_max_longitudinal_performance(node_0, v_0);

    V_dot_max = std::min(V_dot_max, V_dot_x_Vmax);
    //
    const real distance_V_dot_max_1 = this->signed_distance(node_1, this->eval_V_next(node_0, cell, V_dot_max), V_dot_max);
    //
    real VDOTSOL = 0.0;
    if ((distance_V_dot_max_1 <= 0) )
    {
      VDOTSOL = V_dot_max;
    }
    else 
    {
      this->BD.solve(
        [this,&node_0, &node_1, &cell](const real VDOT) -> real {
          return this->signed_distance(node_1, this->eval_V_next(node_0, cell, VDOT), VDOT) + this->solver_p.tolerance; },
        V_dot_0,
        V_dot_max,
        VDOTSOL );
    }
    // check if the solution is valid
    if (!std::isnan(VDOTSOL))
    {
      // update cell
      cell.V_dot = VDOTSOL;
      cell.V_1   = this->eval_V_next(node_0, cell, VDOTSOL);
      // update next v0 with the final velocity of the segment of the maximum reachable velocity
      v_0 = cell.V_1;
      // set type
      cell.m_type = SegmentType::FORWARD;
      // set Vmax of the following segment to be the one of the current cell
      node_1.V_max = std::min(node_1.V_max, cell.V_1);
    }
    else
    {
      // set type
      cell.m_type = SegmentType::FORWARD_NAN;
      this->dump_seg_id.push_back(cell.ID0);
      //set the next vmax
      node_1.V_max = std::min(node_1.V_max, v_0);
      // set v_0 to be the maximum reachable velocity
      v_0 = node_1.V_max;
      cell.V_1 = v_0;
    }
  }
}

// --------------------------------------------------------------------------------------------


void
FBGA_MOTO::BW()
{
  real v_1 = this->Nodes.back().V_max;
  // if last segment is valid forward ten use the final setted velocity
  if(!std::isnan(this->Cells.back().V_dot))
  {
    v_1 = this->Cells.back().V_1;
  }
  // consider each cell but in reverse order
  for (auto i = static_cast<integer>(this->Cells.size() - 1); i >= 0; i--)
  {
    auto & cell = this->Cells[i];
    auto & node_0 = this->Nodes[cell.ID0];
    auto & node_1 = this->Nodes[cell.ID1];
    // get final velocity. If final index final velocity is the minimum of the reached velocity and
    // the final velocity of the segment. Otherwise, it is the initial velocity of the next segment.
    cell.V_1 = v_1;
    // ge the max and min final acceleration 
    const real V_dot_max = this->get_max_longitudinal_performance(node_1, v_1);
    const real V_dot_min = this->get_min_longitudinal_performance(node_1, v_1);

    const real v0_reach_max = std::min(node_0.V_max, this->eval_V_prev(node_1, cell, V_dot_min));
    const real v0_reach_min = std::max(0.0,          this->eval_V_prev(node_1, cell, V_dot_max));
    // 
    // check if cell.V_0 is inside the reachable velocity range
    const bool is_v0_reachable = (cell.V_0 >= v0_reach_min && cell.V_0 <= v0_reach_max);
    // compute the mean acceleration
    const real a_mean = (cell.V_1 * cell.V_1 - cell.V_0 * cell.V_0) / (2.0 * cell.L);
    // const real V_dot_mean = this->eval_V_dot_Vahatx(node_1, cell.V_1, a_mean);
    const real V_dot_mean = this->eval_V_dot_Vahatx(node_0, cell.V_0, a_mean);
    // check if V_dot_mean is inside the V_dot_min - V_dot_max range
    const bool is_V_dot_mean_candidate = ((V_dot_mean >= V_dot_min) && (V_dot_mean <= V_dot_max));
    // check if the final point is still feasible with V_dot_mean
    const real distance_0_mean = this->signed_distance(node_0, this->eval_V_prev(node_0, cell, V_dot_mean), V_dot_mean);
    // if the distance is less than the tolerance, we can use V_dot_mean and a candidate
    // if the segment is a forward pass and the final saved velocity coincides with v_1 ignore and continue
    const bool is_valid_forward = (cell.m_type == SegmentType::FORWARD && (std::abs(this->eval_V_next(node_0, cell, cell.V_dot) - v_1) <= this->solver_p.tolerance));
    const bool is_valid_yellow_flag = (cell.m_type == SegmentType::YELLOWFLAG && (std::abs(this->eval_V_next(node_0, cell, cell.V_dot) - v_1) <= this->solver_p.tolerance));
    if(is_valid_forward || is_valid_yellow_flag)
    {
      v_1 = cell.V_0; // update v_1 to the current V_0
      continue;
    }
    if (is_V_dot_mean_candidate && is_v0_reachable && (distance_0_mean <= (2*this->solver_p.tolerance)) )
    {
      // set the V_dot to the mean value
      cell.V_dot = V_dot_mean;
      //cell.V_0   = this->eval_V_prev(node_0, cell, V_dot_mean);
      cell.m_type = SegmentType::TRANSITION;
      // update v_1
      v_1 = cell.V_0;
      continue; // continue with the next cell
    }
    else
    {
      const real distance_amin = this->signed_distance(node_0, this->eval_V_prev(node_0, cell, V_dot_min), V_dot_min);
      real VDOTSOL = 0.0;
      if (distance_amin <= 0)
      {
        VDOTSOL = V_dot_min;
      }
      else 
      {
        // choose a for the solver (with a value inside the envelope)
        const real V_dot_solver = this->eval_V_dot_Vatildex(node_1, cell.V_1, 0.0);
        this->BD.solve(
          [this,&node_0, &cell](const real VDOT) -> real {
            return this->signed_distance(node_0, this->eval_V_prev(node_0, cell, VDOT), VDOT) + this->solver_p.tolerance; },
          V_dot_min,
          V_dot_solver,
          VDOTSOL );
      }
      if (!std::isnan(VDOTSOL))
      {
        // set v dot and change the v0
        cell.V_dot = VDOTSOL;
        cell.V_0   = this->eval_V_prev(node_0, cell, VDOTSOL);
        cell.m_type = SegmentType::BACKWARD;
        //update v_1
        v_1 = cell.V_0;
      }
      else 
      {
        std::cout << "FBGA_MOTO::BW() >> No solution found for cell " << i << "\n";
        // set type
        cell.V_dot = QUIET_NAN;
        cell.m_type = SegmentType::BACKWARD_NAN;
        this->dump_seg_id.push_back(cell.ID0);
      }
    }
  }
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::compute_total_time() const
{
  real T = 0;
  for(const auto & cell: this->Cells)
  {
    T += static_cast<real>(2)*cell.L / (cell.V_0 + cell.V_1);
  }
  return T;
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::compute_time(real s) const
{
  real T = 0;
  for(const auto & cell: this->Cells)
  {
    if (cell.s_0 > s)
    {
      break;
    }
    if ((s > cell.s_0) && (s < cell.s_1))
    {
      const auto & node = this->Nodes[cell.ID0];
      const real s_loc = s - cell.s_0;
      const real V_1 = this->eval_V(node,cell,s_loc);
      T += static_cast<real>(2) * s_loc / (cell.V_0 + V_1);
    }
    else 
    {
      T += static_cast<real>(2) * cell.L / (cell.V_0 + cell.V_1);
    }
  }
  return T;
}

// --------------------------------------------------------------------------------------------

inline void
FBGA_MOTO::compute_context(const GG::NodeStruct3D& node, real V, real V_dot, Context & ctx) const
{
  // Basic trigonometric and kinematic calculations
    ctx.cos_chi = std::cos(node.chi);
    ctx.sin_chi = std::sin(node.chi);
    ctx.den_common = node.Omega_z * node.n - 1.0;
    ctx.inv_den_common = 1.0 / ctx.den_common;
    ctx.inv_den_common2 = std::pow(ctx.inv_den_common, 2.0);
    ctx.inv_den_common3 = std::pow(ctx.inv_den_common, 3.0);
    
    // Velocity calculations
    ctx.s_dot = V * ctx.cos_chi * (-ctx.inv_den_common);
    ctx.w = node.Omega_x * ctx.s_dot * node.n;
    ctx.chi_dot = this->m_ablation_flags.compute_chi_dot ? node.chi_prime * ctx.s_dot : 0.0;
    ctx.n_dot = V * ctx.sin_chi;
    
    // Angular velocity transformations
    ctx.omega_hat_x = (ctx.sin_chi * node.Omega_y + ctx.cos_chi * node.Omega_x) * ctx.s_dot;
    ctx.omega_hat_y = (ctx.cos_chi * node.Omega_y - ctx.sin_chi * node.Omega_x) * ctx.s_dot;
    ctx.omega_hat_z = node.Omega_z * ctx.s_dot + ctx.chi_dot;
    
    // Second derivative calculations
    // ctx.s_ddotA = V * V * (-node.Omega_z_prime * ctx.cos_chi * node.n + node.Omega_z * ctx.cos_chi * ctx.den_common) * ctx.cos_chi * ctx.inv_den_common3;
    ctx.s_ddotA = - V * ctx.cos_chi * ( node.Omega_z_prime * ctx.s_dot * node.n + node.Omega_z * ctx.sin_chi * V) * ctx.inv_den_common2;
    ctx.s_ddotB = (ctx.chi_dot * ctx.cos_chi * V - ctx.cos_chi * V_dot) * ctx.inv_den_common;
    ctx.s_ddot = (!this->m_ablation_flags.compute_s_ddot) ? 0.0 : ctx.s_ddotA + ctx.s_ddotB;
    ctx.w_dot = node.Omega_x * ctx.s_dot * ctx.n_dot + node.Omega_x * ctx.s_ddot * node.n + node.Omega_x_prime * ctx.s_dot * ctx.s_dot * node.n; // is it supposed to be s_dot squared?
    
    // Acceleration components
    ctx.a_tilde_x = (ctx.omega_hat_y * ctx.w) + V_dot + node.g_x;
    ctx.a_tilde_y = (ctx.omega_hat_z * V) - (ctx.omega_hat_x * ctx.w) + node.g_y;
    ctx.a_tilde_z = ctx.w_dot - (ctx.omega_hat_y * V) + node.g_z;
    
    // G-G diagram parameters
    ctx.a_tilde_y_lim = (node.alpha * this->gggv_moto.a_y_lim(V, ctx.a_tilde_z) - this->m_lat_tol);
    ctx.a_tilde_y_clip = clip(ctx.a_tilde_y, -ctx.a_tilde_y_lim, ctx.a_tilde_y_lim);
    
    ctx.a_tilde_x_max_gg = node.alpha * this->gggv_moto.a_x_push(ctx.a_tilde_y_clip, V, ctx.a_tilde_z);
    ctx.a_tilde_x_min_gg = node.alpha * this->gggv_moto.a_x_pull(ctx.a_tilde_y_clip, V, ctx.a_tilde_z);
    ctx.a_tilde_x_eng = this->gggv_moto.a_x_eng(V);
    ctx.rho_max = 0.0;
    ctx.rho_min = 0.0;
    
    // Compute final bounds
    ctx.a_tilde_x_max = std::min(
        ctx.a_tilde_x_max_gg,
        ctx.a_tilde_x_eng
    );
    
    ctx.a_tilde_x_min = ctx.a_tilde_x_min_gg;
    
    // Handle edge case at lateral limit
    ctx.at_lateral_limit = (std::abs(ctx.a_tilde_y_clip) >= ctx.a_tilde_y_lim);
    if (ctx.at_lateral_limit) {
        ctx.a_tilde_x_max = std::min(0.0, ctx.a_tilde_x_eng);
        ctx.a_tilde_x_min = std::max(0.0, ctx.a_tilde_x_min);
    }
    
    // Alternative check for tolerance-based edge case (used in min/max functions)
    if ((std::abs(ctx.a_tilde_y) - ctx.a_tilde_y_lim) > this->solver_p.tolerance) {
        ctx.a_tilde_x_max = std::min(0.0, ctx.a_tilde_x_eng);
        ctx.a_tilde_x_min = std::max(0.0, ctx.a_tilde_x_min);
    }

}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::signed_distanceNOOPTI(GG::NodeStruct3D const & node, real V, real V_dot)
{
  const real a_tilde_y = this->eval_a_tilde_y(node, V);
  const real a_tilde_lim_y = (this->eval_a_tilde_y_lim(node, V, V_dot) - this->m_lat_tol);
  const real a_tilde_x_max = this->eval_a_tilde_x_max(node, V, V_dot);
  const real a_tilde_x_min = this->eval_a_tilde_x_min(node, V, V_dot);
  const real a_tilde_x = this->eval_a_tilde_x(node, V, V_dot);
  return GG::signed_distance(a_tilde_x, a_tilde_x_min, a_tilde_x_max,
                             a_tilde_y, -a_tilde_lim_y, a_tilde_lim_y);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::signed_distance(GG::NodeStruct3D const & node, real V, real V_dot)
{
  Context ctx;
  this->compute_context(node, V, V_dot, ctx);
  return GG::signed_distance(ctx.a_tilde_x, ctx.a_tilde_x_min, ctx.a_tilde_x_max,
                              ctx.a_tilde_y, -ctx.a_tilde_y_lim, ctx.a_tilde_y_lim);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::signed_distance_YF(GG::NodeStruct3D const & node, real V, real V_dot)
{
  const real a_tilde_y = this->eval_a_tilde_y(node, V);
  const real a_tilde_lim_y = (this->eval_a_tilde_y_lim(node, V, V_dot) - this->m_lat_tol);
  const real a_tilde_x_max = this->eval_a_tilde_x_max(node, V, V_dot);
  const real a_tilde_x_min = this->eval_a_tilde_x_min_YF(node, V, V_dot);
  const real a_tilde_x = this->eval_a_tilde_x(node, V, V_dot);
  return GG::signed_distance(a_tilde_x, a_tilde_x_min, a_tilde_x_max,
                             a_tilde_y, -a_tilde_lim_y, a_tilde_lim_y);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_V_next(GG::NodeStruct3D const & node, GG::CellStruct3D const & cell, real VDOT)
{
  const real a_hat_x = eval_a_hat_x(node, cell.V_0, VDOT);
  return std::sqrt(
    (static_cast<real>(2) * cell.L * a_hat_x) + std::pow(cell.V_0, 2)
  );
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_V_prev(GG::NodeStruct3D const & node, GG::CellStruct3D const & cell, real VDOT)
{
  const real a_hat_x = eval_a_hat_x(node, cell.V_1, VDOT);
  return std::sqrt(
    - static_cast<real>(2) * cell.L * a_hat_x + std::pow(cell.V_1, 2)
  );
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_a_tilde_x(GG::NodeStruct3D const & node, real V, real V_dot)
{
  return eval_a_hat_x(node, V, V_dot) + node.g_x;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_a_hat_x(GG::NodeStruct3D const & node, real V, real V_dot)
{
  const real cos_chi = std::cos(node.chi);
  // const real sin_chi = std::sin(node.chi);
  const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
  // omega_y = (cos(chi) * Omega_y - sin(chi) * Omega_x) * s_dot;
  const real w = node.Omega_x*s_dot*node.n;
  const real omega_hat_y = (std::cos(node.chi) * node.Omega_y - std::sin(node.chi) * node.Omega_x) * s_dot;
  return (omega_hat_y*w) + V_dot;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_V_dot_Vatildex(GG::NodeStruct3D const & node, real V, real atildex)
{
  // omega_x = (sin(chi) * Omega_y + cos(chi) * Omega_x) * s_dot;
  // omega_y = (cos(chi) * Omega_y - sin(chi) * Omega_x) * s_dot;
  // omega_z = Omega_z * s_dot + chi_dot;
  const real cos_chi = std::cos(node.chi);
  // const real sin_chi = std::sin(node.chi);
  const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
  const real w = node.Omega_x*s_dot*node.n;
  const real omega_hat_y = (std::cos(node.chi) * node.Omega_y - std::sin(node.chi) * node.Omega_x) * s_dot;
  return atildex - node.g_x - (omega_hat_y * w);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_V_dot_Vahatx(GG::NodeStruct3D const & node, real V, real ahatx)
{
  // omega_x = (sin(chi) * Omega_y + cos(chi) * Omega_x) * s_dot;
  // omega_y = (cos(chi) * Omega_y - sin(chi) * Omega_x) * s_dot;
  // omega_z = Omega_z * s_dot + chi_dot;
  const real cos_chi = std::cos(node.chi);
  const real sin_chi = std::sin(node.chi);
  const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
  const real w = node.Omega_x*s_dot*node.n;
  const real omega_hat_y = (cos_chi * node.Omega_y - sin_chi * node.Omega_x) * s_dot;
  return ahatx - (omega_hat_y * w);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_a_tilde_y(GG::NodeStruct3D const & node, real V)
{
  return this->eval_a_hat_y(node,V) + node.g_y;
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_a_tilde_z(GG::NodeStruct3D const & node, real V, real V_dot)
{
  return eval_a_hat_z(node,V, V_dot) + node.g_z;
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_a_hat_y(GG::NodeStruct3D const & node, real V)
{
  // omega_x = (sin(chi) * Omega_y + cos(chi) * Omega_x) * s_dot;
  // omega_y = (cos(chi) * Omega_y - sin(chi) * Omega_x) * s_dot;
  // omega_z = Omega_z * s_dot + chi_dot;
  const real cos_chi = std::cos(node.chi);
  const real sin_chi = std::sin(node.chi);
  const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
  const real w = node.Omega_x*s_dot*node.n;
  const real chi_dot = this->m_ablation_flags.compute_chi_dot ? node.chi_prime * s_dot : 0.0; 
  const real omega_hat_x = (sin_chi * node.Omega_y + cos_chi * node.Omega_x) * s_dot;
  const real omega_hat_z = node.Omega_z * s_dot + chi_dot;
  return (omega_hat_z * V) - (omega_hat_x * w);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_a_hat_z(GG::NodeStruct3D const & node, real V, real V_dot)
{
  const real cos_chi = std::cos(node.chi);
  const real sin_chi = std::sin(node.chi);
  const real s_dot = V * cos_chi / (1 - node.n * node.Omega_z);
  const real omega_hat_y = (cos_chi * node.Omega_y - sin_chi * node.Omega_x) * s_dot;
  const real n_dot = V * sin_chi;
  const real s_ddot = this->eval_s_ddot(node,V,V_dot);
  const real w_dot = 
    node.Omega_x * s_dot * n_dot +
    node.Omega_x * s_ddot * node.n +
    node.Omega_x_prime * s_dot * s_dot * node.n; // is it supposed to be s_dot squared?
  return w_dot - (omega_hat_y * V);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_a_tilde_y_lim(GG::NodeStruct3D const & node, real V, real V_dot)
{
  const real a_tilde_z = eval_a_tilde_z(node, V, V_dot);
  return node.alpha * this->gggv_moto.a_y_lim(V,a_tilde_z);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_a_tilde_x_max(GG::NodeStruct3D const & node, real V, real V_dot)
{
  // TODO OPTIMIZE
  const real a_tilde_y = eval_a_tilde_y(node, V);
  const real a_tilde_z = eval_a_tilde_z(node, V, V_dot);
  const real a_tilde_y_lim = (node.alpha * this->gggv_moto.a_y_lim(V,a_tilde_z) - this->m_lat_tol);
  const real a_tilde_y_clip = clip(a_tilde_y, -a_tilde_y_lim, a_tilde_y_lim);

  const real a_tilde_x_max = node.alpha * this->gggv_moto.a_x_push(a_tilde_y_clip, V, a_tilde_z);
  const real a_tilde_x_eng = this->gggv_moto.a_x_eng(V); // no alpha scaling. Engine power do not depend on adherence

  // const real rho_max     = 0.0;

  if (std::abs(std::abs(a_tilde_y) - a_tilde_y_lim) < this->solver_p.tolerance)
  {
    // if a_tilde_y is very close to a_y_lim, we can use the a_x_eng value
    // to avoid numerical issues
    return std::min(0.0, a_tilde_x_eng);
  }
  
  return std::min(
    a_tilde_x_max,
    a_tilde_x_eng
  );
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_a_tilde_x_min(GG::NodeStruct3D const & node, real V, real V_dot)
{
  // TODO OPTIMIZE
  const real a_tilde_y = eval_a_tilde_y(node, V);
  const real a_tilde_z = eval_a_tilde_z(node, V, V_dot);
  const real a_tilde_y_lim = (node.alpha * this->gggv_moto.a_y_lim(V,a_tilde_z) - this->m_lat_tol);
  const real a_tilde_y_clip = clip(a_tilde_y, -a_tilde_y_lim, a_tilde_y_lim);
  const real a_tilde_x_min = node.alpha * this->gggv_moto.a_x_pull(a_tilde_y_clip, V, a_tilde_z);
  // const real rho_min     = 0.0;

  if (std::abs(std::abs(a_tilde_y) - a_tilde_y_lim) < this->solver_p.tolerance)
  {
    // if a_tilde_y is very close to a_y_lim, we can use the a_x_min value
    // to avoid numerical issues
    return std::max(0.0, a_tilde_x_min);
  }

  return a_tilde_x_min;
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_a_tilde_x_min_YF(GG::NodeStruct3D const & node, real V, real V_dot)
{
  // TODO OPTIMIZE
  const real a_tilde_y = eval_a_tilde_y(node, V);
  const real a_tilde_z = eval_a_tilde_z(node, V, V_dot);
  const real a_tilde_y_lim = (node.alpha * this->gggv_moto.a_y_lim(V,a_tilde_z) - this->m_lat_tol);
  const real a_tilde_y_clip = clip(a_tilde_y, -a_tilde_y_lim, a_tilde_y_lim);
  const real a_tilde_x_min = node.alpha * this->gggv_moto.a_x_pull(a_tilde_y_clip, V, a_tilde_z);
  // const real rho_min     = 0.0;

  if (std::abs(std::abs(a_tilde_y) - a_tilde_y_lim) < this->solver_p.tolerance)
  {
    // if a_tilde_y is very close to a_y_lim, we can use the a_x_min value
    // to avoid numerical issues
    return std::max(0.0, a_tilde_x_min);
  }


  return std::max(
    a_tilde_x_min ,
    this->m_yellow_flag_data.a_des_min);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_s_ddot(GG::NodeStruct3D const & node, real V, real V_dot)
{
  if (!this->m_ablation_flags.compute_s_ddot)
  {
    return 0.0; // if s_ddot is not computed, return 0
  }
  const real chi = node.chi;
  const real cos_chi = std::cos(chi);
  const real sin_chi = std::sin(chi);
  // optimized code
  const real den_common = node.Omega_z * node.n - 1.0;
  const real inv_den_common = 1.0 / den_common;
  // const real inv_den_common3 = pow(inv_den_common, 3.0);
  const real inv_den_common2 = inv_den_common * inv_den_common;
  //
  const real s_dot = V * cos_chi * (-inv_den_common);
  //
  const real chi_dot = this->m_ablation_flags.compute_chi_dot ? node.chi_prime * s_dot : 0.0;
  // separate contribution for clarity
  // const real s_ddotA = V * V * (-node.Omega_z_prime * s_dot * cos_chi * node.n + node.Omega_z * cos_chi * den_common) * cos_chi * inv_den_common3;
  const real s_ddotA = - V * cos_chi * ( node.Omega_z_prime * s_dot * node.n + node.Omega_z * V * sin_chi) * inv_den_common2;
  const real s_ddotB = ( chi_dot * cos_chi * V - cos_chi * V_dot )* inv_den_common;
  return s_ddotA + s_ddotB;
}

// --------------------------------------------------------------------------------------------

//  █████╗ ██╗   ██╗██╗  ██╗██╗██╗     ██╗ █████╗ ██████╗ ██╗   ██╗
// ██╔══██╗██║   ██║╚██╗██╔╝██║██║     ██║██╔══██╗██╔══██╗╚██╗ ██╔╝
// ███████║██║   ██║ ╚███╔╝ ██║██║     ██║███████║██████╔╝ ╚████╔╝ 
// ██╔══██║██║   ██║ ██╔██╗ ██║██║     ██║██╔══██║██╔══██╗  ╚██╔╝  
// ██║  ██║╚██████╔╝██╔╝ ██╗██║███████╗██║██║  ██║██║  ██║   ██║   
// ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═╝╚══════╝╚═╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝   


// --------------------------------------------------------------------------------------------

// ██████╗  ██████╗ ███████╗████████╗██████╗ ██████╗  ██████╗  ██████╗███████╗███████╗███████╗██╗███╗   ██╗ ██████╗ 
// ██╔══██╗██╔═══██╗██╔════╝╚══██╔══╝██╔══██╗██╔══██╗██╔═══██╗██╔════╝██╔════╝██╔════╝██╔════╝██║████╗  ██║██╔════╝ 
// ██████╔╝██║   ██║███████╗   ██║   ██████╔╝██████╔╝██║   ██║██║     █████╗  ███████╗███████╗██║██╔██╗ ██║██║  ███╗
// ██╔═══╝ ██║   ██║╚════██║   ██║   ██╔═══╝ ██╔══██╗██║   ██║██║     ██╔══╝  ╚════██║╚════██║██║██║╚██╗██║██║   ██║
// ██║     ╚██████╔╝███████║   ██║   ██║     ██║  ██║╚██████╔╝╚██████╗███████╗███████║███████║██║██║ ╚████║╚██████╔╝
// ╚═╝      ╚═════╝ ╚══════╝   ╚═╝   ╚═╝     ╚═╝  ╚═╝ ╚═════╝  ╚═════╝╚══════╝╚══════╝╚══════╝╚═╝╚═╝  ╚═══╝ ╚═════╝ 
                                                                                                                 
// --------------------------------------------------------------------------------------------


integer FBGA_MOTO::get_cell_idx(const real s) const
{
  // Binary search since cells are sorted by s_0
  integer left = 0;
  integer right = static_cast<integer>(this->Cells.size()) - 1;
  // if s is greater or equal to the last cell's s_1, return the last index
  if (s >= this->Cells[right].s_1) {
    const_cast<FBGA_MOTO*>(this)->m_past_index = right; // Update past_index for next search
    return right;
  }
  // if s is less than the first cell's s_0, return -1 (not found)
  if (s < this->Cells[left].s_0) {
    const_cast<FBGA_MOTO*>(this)->m_past_index = -1; // Update past_index for next search
    return -1;
  }
  
  while (left <= right) {
    const integer mid = left + (right - left) / 2;
    const auto& cell = this->Cells[mid];
    
    if (s >= cell.s_0 && s < cell.s_1) {
      // Update past_index for next search
      const_cast<FBGA_MOTO*>(this)->m_past_index = mid;
      return mid;
    }
    
    if (s < cell.s_0) {
      right = mid - 1;
    } else {
      left = mid + 1;
    }
  }
  
  return -1; // Not found
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_Vmax(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const real s_0 = this->Cells[cell_idx].s_0;
  const real s_1 = this->Cells[cell_idx].s_1;
  const real s_norm = (s - s_0) / (s_1 - s_0);
  const real V_max_0 = this->Nodes[this->Cells[cell_idx].ID0].V_max;
  const real V_max_1 = this->Nodes[this->Cells[cell_idx].ID1].V_max;
  return (V_max_0 * (1 - s_norm)) + (V_max_1 * s_norm);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_Omega_x(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const real s_0 = this->Cells[cell_idx].s_0;
  const real s_1 = this->Cells[cell_idx].s_1;
  const real s_norm = (s - s_0) / (s_1 - s_0);
  const real Omega_x_0 = this->Nodes[this->Cells[cell_idx].ID0].Omega_x;
  const real Omega_x_1 = this->Nodes[this->Cells[cell_idx].ID1].Omega_x;
  return (Omega_x_0 * (1 - s_norm)) + (Omega_x_1 * s_norm);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_Omega_y(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const real s_0 = this->Cells[cell_idx].s_0;
  const real s_1 = this->Cells[cell_idx].s_1;
  const real s_norm = (s - s_0) / (s_1 - s_0);
  const real Omega_y_0 = this->Nodes[this->Cells[cell_idx].ID0].Omega_y;
  const real Omega_y_1 = this->Nodes[this->Cells[cell_idx].ID1].Omega_y;
  return (Omega_y_0 * (1 - s_norm)) + (Omega_y_1 * s_norm);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_Omega_z(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const real s_0 = this->Cells[cell_idx].s_0;
  const real s_1 = this->Cells[cell_idx].s_1;
  const real s_norm = (s - s_0) / (s_1 - s_0);
  const real Omega_z_0 = this->Nodes[this->Cells[cell_idx].ID0].Omega_z;
  const real Omega_z_1 = this->Nodes[this->Cells[cell_idx].ID1].Omega_z;
  return (Omega_z_0 * (1 - s_norm)) + (Omega_z_1 * s_norm);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_V(GG::NodeStruct3D const & node, GG::CellStruct3D const & cell,real s)
{
  const real a_hat_x = eval_a_hat_x(node, cell.V_0, cell.V_dot);
  return std::sqrt(
    (static_cast<real>(2) * s * a_hat_x) + std::pow(cell.V_0, 2)
  );
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_V(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0; 
  return eval_V(node, cell, s_loc);
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_V_dot(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  return cell.V_dot; // V_dot is already computed in the FW method
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_A_hat_x(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0;
  const real V = eval_V(node, cell, s_loc);
  return this->eval_a_hat_x(node, V, cell.V_dot);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_A_hat_y(real s)
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0;
  const real V = eval_V(node, cell, s_loc);
  return this->eval_a_hat_y(node, V);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_A_hat_z(real s)
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0;
  const real V = eval_V(node, cell, s_loc);
  const real V_dot = cell.V_dot; 
  return this->eval_a_hat_z(node, V, V_dot);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_A_tilde_x(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0;
  const real V = eval_V(node, cell, s_loc);
  return this->eval_a_tilde_x(node, V, cell.V_dot);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_A_tilde_y(real s)
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0;
  const real V = eval_V(node, cell, s_loc);
  return this->eval_a_tilde_y(node, V);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_A_tilde_z(real s)
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node = this->Nodes[cell.ID0];
  const real s_loc = s - cell.s_0;
  const real V = eval_V(node, cell, s_loc);
  const real V_dot = cell.V_dot; // V_dot is already computed in the FW method
  return this->eval_a_tilde_z(node, V, V_dot);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_g_x(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node_0 = this->Nodes[cell.ID0];
  const auto & node_1 = this->Nodes[cell.ID1];
  const real s_norm = (s - cell.s_0) / cell.L;
  return (node_0.g_x * (1 - s_norm)) + (node_1.g_x * s_norm);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_g_y(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node_0 = this->Nodes[cell.ID0];
  const auto & node_1 = this->Nodes[cell.ID1];
  const real s_norm = (s - cell.s_0) / cell.L;
  return (node_0.g_y * (1 - s_norm)) + (node_1.g_y * s_norm);
}

// --------------------------------------------------------------------------------------------

real
FBGA_MOTO::eval_g_z(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & cell = this->Cells[cell_idx];
  const auto & node_0 = this->Nodes[cell.ID0];
  const auto & node_1 = this->Nodes[cell.ID1];   
  const real s_norm = (s - cell.s_0) / cell.L;
  return (node_0.g_z * (1 - s_norm)) + (node_1.g_z * s_norm);
}

// --------------------------------------------------------------------------------------------

SegmentType 
FBGA_MOTO::eval_segment_type(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  return this->Cells[cell_idx].m_type;
}

// --------------------------------------------------------------------------------------------

real 
FBGA_MOTO::eval_alpha(real s) const
{
  const integer cell_idx = this->get_cell_idx(s);
  const auto & node_0 = this->Nodes[this->Cells[cell_idx].ID0];
  const auto & node_1 = this->Nodes[this->Cells[cell_idx].ID1];
  // compute the alpha as the maximum of the two nodes
  const real alpha_0 = node_0.alpha;
  const real alpha_1 = node_1.alpha;
  const real s_loc = s - this->Cells[cell_idx].s_0;
  const real s_norm = s_loc / this->Cells[cell_idx].L;
  // linear interpolation
  return (alpha_0 * (1 - s_norm)) + (alpha_1 * s_norm);
}

// --------------------------------------------------------------------------------------------

output_plot_ggv_shell 
FBGA_MOTO::eval_shell_plot(GG::NodeStruct3D const & fictitious_node, const real Va, const real Vb, const real V_dot
) 
{
  output_plot_ggv_shell output_plot;
  // for velocity (20 points) from Va to Vb compute a_lim. then another for -alim to +alim
  constexpr size_type num_points = NUMPTSPLOTSHELL;
  const real delta_v = (Vb - Va) / static_cast<real>(num_points - 1);
  for (size_type i = 0; i < num_points; ++i)
  {
    const real V = Va + (static_cast<real>(i) * delta_v);
    // compute a_tilde_z
    const real a_tilde_z = this->eval_a_tilde_z(fictitious_node, V, V_dot);
    // compute a_tilde_y_lim
    const real a_tilde_y_lim = fictitious_node.alpha *  this->gggv_moto.a_y_lim(V,a_tilde_z);
    const real delta_a_tilde_y = 2.0 * a_tilde_y_lim / (real)(num_points - 1);
    for(size_type j = 0; j < num_points; ++j)
    {
      const real a_tilde_y      = -a_tilde_y_lim + (static_cast<real>(j) * delta_a_tilde_y);
      const real a_tilde_y_clip = clip(a_tilde_y, -a_tilde_y_lim, a_tilde_y_lim);
      // compute a_tilde_x_max
      const real a_tilde_x_max = this->gggv_moto.a_x_push(a_tilde_y_clip, V, a_tilde_z);
      // compute a_tilde_x_min
      const real a_tilde_x_min = this->gggv_moto.a_x_pull(a_tilde_y_clip, V, a_tilde_z);
      // store the point
      output_plot.a_tilde_x_max.push_back(static_cast<float>(a_tilde_x_max));
      output_plot.a_tilde_x_min.push_back(static_cast<float>(a_tilde_x_min));
      output_plot.a_tilde_y.push_back(static_cast<float>(a_tilde_y));
      output_plot.v.push_back(static_cast<float>(V));
    }
  }
  return output_plot;
}

// --------------------------------------------------------------------------------------------

output_plot_ggv_shell 
FBGA_MOTO::eval_shell_plotpy(GG::input_plot_ggv_shell const & input_plot) 
{
  output_plot_ggv_shell output_plot;
  constexpr size_type num_points = NUMPTSPLOTSHELL;
  for (auto V : input_plot.v)
  {
    const real a_tilde_y_lim = input_plot.alpha *  this->gggv_moto.a_y_lim(V,input_plot.az) * 0.99999;
    const real delta_a_tilde_y = 2.0 * a_tilde_y_lim / static_cast<real>(num_points - 1);
    for(size_type j = 0; j < num_points; ++j)
    {
      const real a_tilde_y      = -a_tilde_y_lim + (static_cast<real>(j) * delta_a_tilde_y);
      const real a_tilde_y_clip = clip(a_tilde_y, -a_tilde_y_lim, a_tilde_y_lim);
      // const real rho_max            = 0.0;
      // const real rho_min            = 0.0;
      // compute a_tilde_x_max
      const real a_tilde_x_max = std::min(
        input_plot.alpha * this->gggv_moto.a_x_push(a_tilde_y_clip, V, input_plot.az), 
        this->gggv_moto.a_x_eng(V)
      );
      // compute a_tilde_x_min
      const real a_tilde_x_min = input_plot.alpha * this->gggv_moto.a_x_pull(a_tilde_y_clip, V, input_plot.az);
      // store the point
      output_plot.a_tilde_x_max.push_back(static_cast<float>(a_tilde_x_max));
      output_plot.a_tilde_x_min.push_back(static_cast<float>(a_tilde_x_min));
      output_plot.a_tilde_y.push_back(static_cast<float>(a_tilde_y));
      output_plot.v.push_back(static_cast<float>(V));
    }
  }
  return output_plot;
}

// --------------------------------------------------------------------------------------------

void 
FBGA_MOTO::eval_contraint_satisfaction()
{
  // This function is used to evaluate the constraint satisfaction for the FBGA_INDY method.
  // It computes the necessary parameters and updates the context accordingly.
  Context ctx0;
  Context ctx1;
  real cs0 = 0.0;
  real cs1 = 0.0;
  real tmp1 = 0.0;
  real tmp2 = 0.0;
  real tmp3 = 0.0;
  real tmp4 = 0.0;
  //
  std::cout << "Evaluating constraint satisfaction for FBGA_INDY..." << std::endl;
  for (const auto& cell : this->Cells) 
  {
    const auto& node0 = this->Nodes[cell.ID0];
    const auto& node1 = this->Nodes[cell.ID1];
    const real V0 = cell.V_0; // V_0 is already computed in the FW method
    const real V1 = cell.V_1; // V_1 is already computed in
    const real V_dot = cell.V_dot; // V_dot is already computed in the FW method

    // Compute context for the current node
    this->compute_context(node0, V0, V_dot, ctx0);
    this->compute_context(node1, V1, V_dot, ctx1);

    if (ctx0.a_tilde_x >= 0)
    {
      tmp1 = std::pow((ctx0.a_tilde_x / ctx0.a_tilde_x_max_gg), ctx0.rho_max);
      tmp2 = std::pow((std::abs(ctx0.a_tilde_y_clip) / ctx0.a_tilde_y_lim), ctx0.rho_max);
      cs0 = -1 + tmp1 + tmp2; 
    }
    else
    {
      tmp1 = std::pow((ctx0.a_tilde_x / ctx0.a_tilde_x_min_gg), ctx0.rho_min);
      tmp2 = std::pow((std::abs(ctx0.a_tilde_y_clip) / ctx0.a_tilde_y_lim), ctx0.rho_min);
      cs0 = -1 + tmp1 + tmp2;
    }

    if( ctx1.a_tilde_x >= 0)
    {
      tmp3 = std::pow((ctx1.a_tilde_x / ctx1.a_tilde_x_max_gg), ctx1.rho_max);
      tmp4 = std::pow((std::abs(ctx1.a_tilde_y_clip) / ctx1.a_tilde_y_lim), ctx1.rho_max);
      cs1 = -1 + tmp3 + tmp4;
    }
    else
    {
      tmp3 = std::pow((ctx1.a_tilde_x / ctx1.a_tilde_x_min_gg), ctx1.rho_min);
      tmp4 = std::pow((std::abs(ctx1.a_tilde_y_clip) / ctx1.a_tilde_y_lim), ctx1.rho_min);
      cs1 = -1 + tmp3 + tmp4;
    }
    std::cout << "Cell: " << cell.ID0 << " - CS0: " << cs0 << ", CS1: " << cs1 << std::endl;
    // std::cout << "\ttmp1: " << tmp1 << ", tmp2: " << tmp2 << ", tmp3: " << tmp3 << ", tmp4: " << tmp4 << std::endl;
    // std::cout << "\ta_tilde_x: " << ctx0.a_tilde_x << ", a_tilde_y: " << ctx0.a_tilde_y << ", a_tilde_z: " << ctx0.a_tilde_z << std::endl;
    // std::cout << "\ta_tilde_x_max: " << ctx0.a_tilde_x_max << ", a_tilde_x_min: " << ctx0.a_tilde_x_min << std::endl;
    // std::cout << "\ta_tilde_y_lim: " << ctx0.a_tilde_y_lim << ", a_tilde_y_clip: " << ctx0.a_tilde_y_clip << std::endl;
  

    
  }
}


