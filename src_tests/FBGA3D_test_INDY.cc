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

#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <thread>
#include <vector>

#include "FBGA3D/FBGA_INDY.hxx"
#include "FBGA3D/gg_utils.hxx"
#include "FBGA3D/types.hxx"


// rapidcsv
#include "rapidcsv.h"


#include "json.hpp"

#include "cxxopts.hpp"

#include "npy.hpp"

// chrono

using Colour = struct Colour {
  uint8_t r{0};
  uint8_t g{0};
  uint8_t b{0};
};

std::map<std::string, Colour> colours = {
    {"white", {255, 255, 255}},    {"black", {0, 0, 0}},
    {"red", {255, 0, 0}},          {"green", {0, 255, 0}},
    {"blue", {0, 0, 255}},         {"yellow", {255, 255, 0}},
    {"cyan", {0, 255, 255}},       {"magenta", {255, 0, 255}},
    {"grey", {128, 128, 128}},     {"orange", {255, 165, 0}},
    {"purple", {128, 0, 128}},     {"brown", {165, 42, 42}},
    {"pink", {255, 192, 203}},     {"turquoise", {64, 224, 208}},
    {"gold", {255, 215, 0}},       {"silver", {192, 192, 192}},
    {"lime", {0, 255, 0}},         {"maroon", {128, 0, 0}},
    {"olive", {128, 128, 0}},      {"teal", {0, 128, 128}},
    {"navy", {0, 0, 128}},         {"beige", {245, 245, 220}},
    {"lavender", {230, 230, 250}}, {"salmon", {250, 128, 114}},
    {"ivory", {255, 255, 240}},    {"khaki", {240, 230, 140}},
    {"indigo", {75, 0, 130}},      {"azure", {240, 255, 255}},
    {"mint", {189, 252, 201}},     {"apricot", {251, 206, 177}},
    {"coral", {255, 127, 80}},     {"crimson", {220, 20, 60}},
    {"fuchsia", {255, 0, 255}},    {"orchid", {218, 112, 214}},
    {"plum", {221, 160, 221}},     {"snow", {255, 250, 250}}};

int main(int argc, char *argv[])
{
  // Parse command line arguments
  cxxopts::Options options("GIGI");
  options.add_options()
    ("h,help", "Print help")
    ("c,circuit", "Circuit name", cxxopts::value<std::string>()->default_value("Yas_Marina_raceline.csv"))
    ("f,folder", "Folder path of the circuit data", cxxopts::value<std::string>()->default_value("./data/INDY/"))
    ("y,yellow", "Yellow flag active", cxxopts::value<bool>()->default_value("false"))
    ("s,speed", "Target Speed yellow flag", cxxopts::value<GG::real>()->default_value("100"))
    ("a,acceleration", "Target acceleration yellow flag", cxxopts::value<GG::real>()->default_value("-5"))
    ("d,dimension", "2D or 3D (default false = 2D)", cxxopts::value<bool>()->default_value("false"))
    ("o,output", "Output file write", cxxopts::value<bool>()->default_value("false"))
    ("n,note", "Added note to output name", cxxopts::value<std::string>()->default_value(""))
    ("b,bound", "restrict and change the engine envelope and scaling factor", cxxopts::value<bool>()->default_value("false"))
    ("m,mesh", "Mesh flag", cxxopts::value<bool>()->default_value("false"));
  
  auto result = options.parse(argc, argv); 
  // Print help
  if (result.count("help"))
  {
    std::cout << options.help() << std::endl;
    exit(0);
  }

  const bool yellow_flag_active = result["yellow"].as<bool>();
  const bool is_2d = result["dimension"].as<bool>();
  const std::string circuit_name = result["circuit"].as<std::string>();
  const std::string folder_path = result["folder"].as<std::string>();
  const bool output_write = result["output"].as<bool>();

  const GG::real v_des_yf = result["speed"].as<GG::real>();
  const GG::real a_des_yf = result["acceleration"].as<GG::real>();

  const std::string note = result["note"].as<std::string>();

  const bool bound = result["bound"].as<bool>();
  const bool mesh_flag = result["mesh"].as<bool>();

  // 

  std::cout << "FBGA3D test started\n";
  std::cout << "Example class FBGA3D\n";

  std::cout << "Is 2D: " << (is_2d ? "Yes" : "No") << "\n";
  std::cout << "Circuit name: " << circuit_name << "\n";
  std::cout << "Folder path: " << folder_path << "\n";
  std::cout << "Yellow flag active: " << (yellow_flag_active ? "Yes" : "No") << "\n";
  std::cout << "Target Speed yellow flag: " << v_des_yf << " m/s\n";
  std::cout << "Target acceleration yellow flag: " << a_des_yf << " m/s²\n";
  std::cout << "Output write: " << (output_write ? "Yes" : "No") << "\n";
  std::cout << "Note: " << note << "\n";

  // ██████╗  ██████╗  █████╗ ██████╗
  // ██╔══██╗██╔═══██╗██╔══██╗██╔══██╗
  // ██████╔╝██║   ██║███████║██║  ██║
  // ██╔══██╗██║   ██║██╔══██║██║  ██║
  // ██║  ██║╚██████╔╝██║  ██║██████╔╝
  // ╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═╝╚═════╝

  const std::string road_file_path = folder_path + circuit_name;


  //
  rapidcsv::Document doc(
    road_file_path, 
    rapidcsv::LabelParams(0, -1), // Ensure the first non-comment line is treated as headers
    rapidcsv::SeparatorParams(',',true),
    rapidcsv::ConverterParams(),
    rapidcsv::LineReaderParams(true, '#')); // Skip only comment lines starting with '#'
  // get column names from csv file
  std::vector<std::string> column_names = {
    "n_rl_m", 
    "chi_rl_rad", 
    "s_ref_rl_m", 
    "theta_ref_rl_rad", 
    "mu_ref_rl_rad", 
    "phi_ref_rl_rad", 
    "dtheta_ref_rl_radpm", 
    "dmu_ref_rl_radpm", 
    "dphi_ref_rl_radpm"
  };

  //
  std::vector<std::vector<GG::real>> YVEC_SPLNE;
  YVEC_SPLNE.reserve(column_names.size());
  for (const auto& col_name : column_names) {
    std::vector<GG::real> TMP = doc.GetColumn<GG::real>(col_name);
    if (TMP.empty()) {
      std::cout << "Column " << col_name << " is empty.\n";
    } else {
      auto it = std::find_if(TMP.begin(), TMP.end(), [](GG::real x) { return std::isnan(x); });
      //resize all vectors to avoid all nans
      std::size_t index = std::distance(TMP.begin(), it);
      TMP.resize(index);
      YVEC_SPLNE.push_back(TMP);
    }
  }
  //
  std::vector<GG::real> Sref = doc.GetColumn<GG::real>("s_ref_rl_m");
  auto it = std::find_if(Sref.begin(), Sref.end(), [](GG::real x) { return std::isnan(x); });
  //resize all vectors to avoid all nans
  std::size_t index = std::distance(Sref.begin(), it);
  Sref.resize(index);
  for (auto& vec : YVEC_SPLNE) {
    vec.resize(index);
  }
  //
  GG::LinearInterpolatorSet traj_spline(
    Sref, 
    YVEC_SPLNE, // X
    column_names);
  ///
  GG::real total_L = Sref.back(); 
  //
  // ███████╗██╗    ██╗██████╗ ██╗    ██╗
  // ██╔════╝██║    ██║██╔══██╗██║    ██║
  // █████╗  ██║ █╗ ██║██████╔╝██║ █╗ ██║
  // ██╔══╝  ██║███╗██║██╔══██╗██║███╗██║
  // ██║     ╚███╔███╔╝██████╔╝╚███╔███╔╝
  // ╚═╝      ╚══╝╚══╝ ╚═════╝  ╚══╝╚══╝
  //

  const std::vector<GG::real> g_data     =std::move((npy::read_npy<GG::real>(folder_path+"g_list.npy")).data);
  const std::vector<GG::real> v_data     =std::move((npy::read_npy<GG::real>(folder_path+"v_list.npy")).data);
  const std::vector<GG::real> ay_max_data=std::move((npy::read_npy<GG::real>(folder_path+"ay_max.npy")).data);
  const std::vector<GG::real> ax_min_data=std::move((npy::read_npy<GG::real>(folder_path+"ax_min.npy")).data);
  const std::vector<GG::real> ax_max_data=std::move((npy::read_npy<GG::real>(folder_path+"ax_max.npy")).data);


  const std::vector<GG::real> v_eng_data = {0.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0};
  std::vector<GG::real> ax_eng_data;
  GG::scaling_gggv_factors scaling_factors;
  if(bound)
  {
    ax_eng_data = {10.0, 10.0,  9.5,  9.0, 7.4, 4.9, 2.9, 0.6, 0.1, 0.0};
    scaling_factors = {
      1.2, // ax_max_scale
      1.35, // ax_min_scale
      1.1, // ay_scale
      1.55, // gg_exponent_ax_pos
      1.80  // gg_exponent_ax_neg
    };
  }
  else 
  {
    ax_eng_data = {5.0, 10.0, 13.0, 10.0, 7.0, 4.6, 2.5, 0.3, 0.0, 0.0};
    scaling_factors = {
      1.0, // ax_max_scale
      1.0, // ax_min_scale
      1.0, // ay_scale
      1.3, // gg_exponent_ax_pos
      1.3  // gg_exponent_ax_neg
    };
  }

  GG::spline_data_collection spline_data {
    v_data,
    g_data,
    ay_max_data,
    ax_min_data,
    ax_max_data,
    scaling_factors,
    v_eng_data,
    ax_eng_data
  };



  GG::FBGA_INDY fbga_indy(spline_data);



  //
  GG::integer numpts_gg = Sref.size();
  std::vector<GG::real> S_gg(numpts_gg);
  std::vector<GG::real> THETA_gg(numpts_gg), PHI_GG(numpts_gg), MU_GG(numpts_gg);
  std::vector<GG::real> dTHETA_gg(numpts_gg), dPHI_GG(numpts_gg), dMU_GG(numpts_gg);
  std::vector<GG::real> CHI_GG(numpts_gg), N_GG(numpts_gg);
  std::vector<GG::real> ALPHAS(numpts_gg); // Adherence coefficients, set to 1.0 for now
  //
  bool activate_alpha = false;
  //
  for (GG::integer i = 0; i < numpts_gg; ++i) {
    GG::real s = Sref[i];
    if(mesh_flag)
    {
      S_gg[i]      = s;
      THETA_gg[i]  = YVEC_SPLNE[doc.GetColumnIdx("theta_ref_rl_rad")][i];
      PHI_GG[i]    = YVEC_SPLNE[doc.GetColumnIdx("phi_ref_rl_rad")][i];
      MU_GG[i]     = YVEC_SPLNE[doc.GetColumnIdx("mu_ref_rl_rad")][i];
      dTHETA_gg[i] = YVEC_SPLNE[doc.GetColumnIdx("dtheta_ref_rl_radpm")][i];
      dPHI_GG[i]   = YVEC_SPLNE[doc.GetColumnIdx("dphi_ref_rl_radpm")][i];
      dMU_GG[i]    = YVEC_SPLNE[doc.GetColumnIdx("dmu_ref_rl_radpm")][i];
      CHI_GG[i]    = YVEC_SPLNE[doc.GetColumnIdx("chi_rl_rad")][i];
      N_GG[i]      = YVEC_SPLNE[doc.GetColumnIdx("n_rl_m")][i];
      ALPHAS[i]    = 1.0;
    }
    else
    {
      S_gg[i]      = s;
      THETA_gg[i]  = traj_spline.eval("theta_ref_rl_rad", s);
      PHI_GG[i]    = traj_spline.eval("phi_ref_rl_rad", s);
      MU_GG[i]     = traj_spline.eval("mu_ref_rl_rad", s);
      dTHETA_gg[i] = traj_spline.eval("dtheta_ref_rl_radpm", s);
      dPHI_GG[i]   = traj_spline.eval("dphi_ref_rl_radpm", s);
      dMU_GG[i]    = traj_spline.eval("dmu_ref_rl_radpm", s);
      CHI_GG[i]    = traj_spline.eval("chi_rl_rad", s);
      N_GG[i]      = traj_spline.eval("n_rl_m", s);
      ALPHAS[i]    = 1.0;
    }
    if(is_2d)
    {
      PHI_GG[i] = 0.0; // Set phi to 0 for 2D case
      MU_GG[i] = 0.0; // Set mu to 0 for
      dPHI_GG[i] = 0.0; // Set dphi to 0 for 2D case
      dMU_GG[i] = 0.0; // Set dmu to
    }
  }
  //
  GG::trajectory_offset_and_angles_container TOA{
      {N_GG, CHI_GG},
      { MU_GG, PHI_GG, THETA_gg,
        dMU_GG, dPHI_GG, dTHETA_gg,
        S_gg},
      {ALPHAS} // Adherence coefficients
    };
  //
  GG::real v_initial = 20.0;
  GG::yellow_flag_data yellow_flag_data{v_des_yf, a_des_yf, yellow_flag_active};

  //  
  auto start = std::chrono::steady_clock::now();
  GG::real T = fbga_indy.compute(TOA, v_initial, yellow_flag_data);
  auto end = std::chrono::steady_clock::now();
  //
  // auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
  auto elapsed_microsec = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << "======================== FBGA STATS ========================\n";
  std::cout << "Execution time: " << (double) elapsed_microsec.count()/1000.0 << " ms\n";
  std::cout << "Total lap time: " << T << " s\n";
  std::cout << "Num segments:   " << numpts_gg << "\n";
  std::cout << "Avg Tcpu/N seg: " << (double) elapsed_microsec.count() / (double) numpts_gg << " microseconds\n";
  std::cout << "=============================================================\n";

  // fbga_indy.eval_contraint_satisfaction();

  if(output_write)
  {
    // Write output to CSV file
    // Remove extension from circuit_name for output file
    std::string base_circuit_name = circuit_name;
    size_t last_dot = base_circuit_name.find_last_of('.');
    if (last_dot != std::string::npos) {
      base_circuit_name = base_circuit_name.substr(0, last_dot);
    }
    std::string output_file_path = folder_path + "FBGA_" + base_circuit_name + "_output_" + note + ".csv";
    std::cout << "Writing output to " << output_file_path << "\n";
    std::ofstream output_file(output_file_path);
    if (output_file.is_open()) {
      std::cout << "Output file opened successfully.\n";
      // preamble
      output_file << "# FBGA output file\n"
        << "# Circuit: " << circuit_name << "\n"
        << "# Folder: " << folder_path << "\n"
        << "# Is 2D: " << (is_2d ? "Yes" : "No") << "\n"
        << "# Yellow flag active: " << (yellow_flag_active ? "Yes" : "No") << "\n"
        << "# Target Speed yellow flag: " << v_des_yf << " m/s\n"
        << "# Target Acceleration yellow flag: " << a_des_yf << " m/s^2\n"

        << "# ======================== FBGA STATS ========================\n"
        << "# Total lap time FBGA: " << T << " s\n"
        << "# Execution time: " << (double) elapsed_microsec.count()/1000.0 << " ms\n"
        << "# Num segments: " << numpts_gg << "\n"
        << "# Avg Tcpu/N seg: " << (double) elapsed_microsec.count() / (double) numpts_gg << " microseconds\n"
        << "# =============================================================\n";

      output_file << "# Solution data:\n";

      output_file << "n_rl_m,chi_rl_rad,s_ref_rl_m," 
      << "theta_ref_rl_rad, mu_ref_rl_rad, phi_ref_rl_rad, dtheta_ref_rl_radpm, dmu_ref_rl_radpm, dphi_ref_rl_radpm,"
      << "alpha_lu,s_fb_m,v_fb_mps,axhat_fb_mps2,ayhat_fb_mps2,azhat_fb_mps2,axtilde_fb_mps2,aytilde_fb_mps2,aztilde_fb_mps2,v_dot_fb_mps2,v_max_fb_mps,g_x_fb_mps2,g_y_fb_mps2,g_z_fb_mps2,Omega_x_fb_radpm,Omega_y_fb_radpm,Omega_z_fb_radpm\n";


      output_file << std::fixed << std::setprecision(12);
      for (GG::integer i = 0; i < numpts_gg; ++i) {
        output_file 
          << traj_spline.eval("n_rl_m", S_gg[i]) << ","
          << traj_spline.eval("chi_rl_rad", S_gg[i]) << ","
          << traj_spline.eval("s_ref_rl_m", S_gg[i]) << ","
          << traj_spline.eval("theta_ref_rl_rad", S_gg[i]) << ","
          << traj_spline.eval("mu_ref_rl_rad", S_gg[i]) << ","
          << traj_spline.eval("phi_ref_rl_rad", S_gg[i]) << ","
          << traj_spline.eval("dtheta_ref_rl_radpm", S_gg[i]) << ","
          << traj_spline.eval("dmu_ref_rl_radpm", S_gg[i]) << ","
          << traj_spline.eval("dphi_ref_rl_radpm", S_gg[i]) << ","
          << fbga_indy.eval_alpha(S_gg[i]) << ","
          << S_gg[i] << ","
          << fbga_indy.eval_V(S_gg[i]) << ","
          << fbga_indy.eval_A_hat_x(S_gg[i]) << ","
          << fbga_indy.eval_A_hat_y(S_gg[i]) << ","
          << fbga_indy.eval_A_hat_z(S_gg[i]) << ","
          << fbga_indy.eval_A_tilde_x(S_gg[i]) << ","
          << fbga_indy.eval_A_tilde_y(S_gg[i]) << ","
          << fbga_indy.eval_A_tilde_z(S_gg[i]) << ","
          << fbga_indy.eval_V_dot(S_gg[i]) << ","
          << fbga_indy.eval_Vmax(S_gg[i]) << ","
      << fbga_indy.eval_g_x(S_gg[i]) << ","
          << fbga_indy.eval_g_y(S_gg[i]) << ","
          << fbga_indy.eval_g_z(S_gg[i]) << ","
          << fbga_indy.eval_Omega_x(S_gg[i]) << ","
          << fbga_indy.eval_Omega_y(S_gg[i]) << ","
          << fbga_indy.eval_Omega_z(S_gg[i]) << "\n";    
      }
    }
    output_file.close();
  }

  


  return 0;
}
