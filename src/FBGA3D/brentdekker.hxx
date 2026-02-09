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

#ifndef BRENTEKKER_HXX
#define BRENTEKKER_HXX

#include "types.hxx"
#include <vector>
#include <string>
#include <functional>


namespace GG
{

  class brentdekker
  {
  private:
    real tol = STD_TOL;
    integer max_iter = STD_MAX_ITER;
    std::string verbose = "iter";
    integer iter = 0;
    bool flag = false;
    real x_last = QUIET_NAN;
    std::vector<real> x_iter;

  public: 
    brentdekker() = default;
    explicit brentdekker(real tol);
    brentdekker(real tol, integer max_iter);
    brentdekker(real tol, integer max_iter, const std::string &verbose);

    void set_tolerance(const real new_tol){ tol = new_tol; }
    void set_max_iter(const integer new_max_iter) { max_iter = new_max_iter; }
    void set_verbose(const std::string &new_verbose) { verbose = new_verbose; }
    void setup(real tol, integer max_iter, const std::string &verbose);
    bool solve(const std::function<real(real)> &f, real a, real b, real & x0);
 }; 

}

#endif
