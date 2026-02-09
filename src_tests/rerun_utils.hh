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
#ifndef INCLUDE_RERUN_UTILS
#define INCLUDE_RERUN_UTILS

#define RERUN_OFF 2
#define RERUN_MINIMAL 1
#define RERUN_ALL 0
#define RERUN_LOGGING RERUN_ALL

#if RERUN_LOGGING < RERUN_OFF
#include "rerun.hpp"
#include "rerun/recording_stream.hpp"
#endif

#if RERUN_LOGGING < RERUN_OFF
#include "rerun/components/scalar.hpp"
#endif
#if RERUN_LOGGING < RERUN_MINIMAL
#include "rerun/archetypes/line_strips3d.hpp"
#include "rerun/collection.hpp"
#include "rerun/components/transform_mat3x3.hpp"
// #include "transforms.hpp"
#endif

#include "FBGA3D/types.hxx"

constexpr float AXIS_RADIUS = 0.05F;

inline void axis_2D(GG::real x_min, GG::real x_max, GG::real y_min, GG::real y_max, const std::string & name)
{
  rerun::RecordingStream &recording = rerun::RecordingStream::current();
  std::vector<std::vector<std::array<float, 2>>> points = {
    {{0.F, (GG::floating) -y_min}, {0.F, (GG::floating) -y_max}},
    {{(GG::floating) x_min, 0.F}, {(GG::floating) x_max, 0.F}}
  };
  recording.log_static(name, 
      rerun::LineStrips2D(points)
      .with_radii({AXIS_RADIUS})
  );
}

#endif