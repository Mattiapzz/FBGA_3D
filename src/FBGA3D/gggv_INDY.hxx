/*
(***********************************************************************************)
(*                                                                                 *)
(* The FBGA3D project                                                               *)
(*                                                                                 *)
(* Copyright (c) 2025, Mattia Piazza                                               *)
(*                                                                                 *)
(* The FBGA3D project and its components are supplied under the terms               *)
(* of the open source BSD 2-Clause License. The contents of the FBGA3D              *)
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
#ifndef GGGV_INDY_HXX
#define GGGV_INDY_HXX

#include "types.hxx"
#include "FBGA3D/gg_utils.hxx"

namespace GG {

  class gggv_INDY { 
    private:

    real m_alpha = 1.0;
    // real m_g = GRAVITY;

    GG::BilinearInterpolator m_ay_max_bilinear;
    GG::BilinearInterpolator m_ax_min_bilinear;
    GG::BilinearInterpolator m_ax_max_bilinear;

    GG::LinearInterpolator m_ax_eng_interpolator;

    scaling_gggv_factors m_scaling_factors;

    public:
    gggv_INDY();
    explicit gggv_INDY(scaling_gggv_factors const & scaling_factors);
    explicit gggv_INDY(spline_data_collection const & spline_data);
    void set_scaling_factors(scaling_gggv_factors const & scaling_factors);
    void set_engine_max(engine_max_ggv const & engine_max_ggv);
    // Pure virtual method to be overridden by derived classes
    [[nodiscard]] real a_x_push( real ay_tilde, real V, real az_tilde) const;
    [[nodiscard]] real a_x_pull( real ay_tilde, real V, real az_tilde) const;
    
    //
    // private:
    void setup_std();
    [[nodiscard]] real rho(real V, real az_tilde) const;
    [[nodiscard]] real rho_max(real V, real az_tilde) const;
    [[nodiscard]] real rho_min(real V, real az_tilde) const;

    [[nodiscard]] real a_x_eng(real V) const;
    [[nodiscard]] real a_x_max(real V, real az_tilde) const;
    [[nodiscard]] real a_x_min(real V, real az_tilde) const;
    [[nodiscard]] real a_y_lim(real V, real az_tilde) const;

  };
  
}
#endif // GGGV_INDY_HXX