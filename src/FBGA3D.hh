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

///
/// file: FBGA3D.hh
///

#pragma once

#ifndef INCLUDE_FBGA3D
#define INCLUDE_FBGA3D


// Print FBGA3D errors
#ifndef FBGA3D_ERROR
#define FBGA3D_ERROR(MSG)                                                                           \
  {                                                                                                \
    throw std::runtime_error(std::to_string(MSG));                                                 \
  }
#endif

// Check for FBGA3D errors
#ifndef FBGA3D_ASSERT
#define FBGA3D_ASSERT(COND, MSG)                                                                    \
  if (!(COND))                                                                                     \
  FBGA3D_ERROR(MSG)
#endif


#endif

///
/// eof: FBGA3D.hh
///
