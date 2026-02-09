#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <rerun/archetypes/points3d.hpp>
#include <rerun/archetypes/scalar.hpp>
#include <thread>
#include <vector>

#include "FBGA3D/FBGA_INDY.hxx"
#include "FBGA3D/gg_utils.hxx"
#include "FBGA3D/types.hxx"

#include "rerun_utils.hh"

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
    ("r,rerun", "rerun active", cxxopts::value<bool>()->default_value("false"))
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

  const bool rerun_active = result["rerun"].as<bool>();
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
  std::cout << "Rerun active: " << (rerun_active ? "Yes" : "No") << "\n";
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
    "x_rl_m", 
    "y_rl_m", 
    "z_rl_m", 
    "v_rl_mps", 
    "n_rl_m", 
    "chi_rl_rad", 
    "ax_rl_mps2", 
    "ay_rl_mps2", 
    "jx_rl_mps3", 
    "jy_rl_mps3", 
    "tire_util_rl", 
    "s_ref_rl_m", 
    "x_ref_rl_m", 
    "y_ref_rl_m" , 
    "z_ref_rl_m", 
    "theta_ref_rl_rad", 
    "mu_ref_rl_rad", 
    "phi_ref_rl_rad", 
    "dtheta_ref_rl_radpm", 
    "dmu_ref_rl_radpm", 
    "dphi_ref_rl_radpm", 
    "w_tr_right_ref_rl_m", 
    "w_tr_left_ref_rl_m", 
    "omega_x_ref_rl_radpm", 
    "omega_y_ref_rl_radpm", 
    "omega_z_ref_rl_radpm"
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
  GG::real v_initial = traj_spline.eval("v_rl_mps", 0.0);
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

      output_file 
        << "x_rl_m,y_rl_m,z_rl_m,v_rl_mps,n_rl_m,chi_rl_rad,ax_rl_mps2,ay_rl_mps2,s_ref_rl_m,x_ref_rl_m,y_ref_rl_m,z_ref_rl_m,omega_x_ref_rl_radpm,omega_y_ref_rl_radpm,omega_z_ref_rl_radpm,"
      
        << "alpha_lu,s_lu_m,v_lu_mps,axhat_lu_mps2,ayhat_lu_mps2,azhat_lu_mps2,axtilde_lu_mps2,aytilde_lu_mps2,aztilde_lu_mps2,v_dot_lu_mps2,v_max_lu_mps,g_x_lu_mps2,g_y_lu_mps2,g_z_lu_mps2,Omega_x_lu_radpm,Omega_y_lu_radpm,Omega_z_lu_radpm\n";


      output_file << std::fixed << std::setprecision(12);
      for (GG::integer i = 0; i < numpts_gg; ++i) {
        output_file 
          << traj_spline.eval("x_rl_m", S_gg[i]) << ","
          << traj_spline.eval("y_rl_m", S_gg[i]) << ","
          << traj_spline.eval("z_rl_m", S_gg[i]) << ","
          << traj_spline.eval("v_rl_mps", S_gg[i]) << ","
          << traj_spline.eval("n_rl_m", S_gg[i]) << ","
          << traj_spline.eval("chi_rl_rad", S_gg[i]) << ","
          << traj_spline.eval("ax_rl_mps2", S_gg[i]) << ","
          << traj_spline.eval("ay_rl_mps2", S_gg[i]) << ","
          << traj_spline.eval("s_ref_rl_m", S_gg[i]) << ","
          << traj_spline.eval("x_ref_rl_m", S_gg[i]) << ","
          << traj_spline.eval("y_ref_rl_m", S_gg[i]) << ","
          << traj_spline.eval("z_ref_rl_m", S_gg[i]) << ","
          << traj_spline.eval("omega_x_ref_rl_radpm", S_gg[i]) << ","
          << traj_spline.eval("omega_y_ref_rl_radpm", S_gg[i]) << ","
          << traj_spline.eval("omega_z_ref_rl_radpm", S_gg[i]) << ","
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

  //
  // ██████╗ ███████╗██████╗ ██╗   ██╗███╗   ██╗
  // ██╔══██╗██╔════╝██╔══██╗██║   ██║████╗  ██║
  // ██████╔╝█████╗  ██████╔╝██║   ██║██╔██╗ ██║
  // ██╔══██╗██╔══╝  ██╔══██╗██║   ██║██║╚██╗██║
  // ██║  ██║███████╗██║  ██║╚██████╔╝██║ ╚████║
  // ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝ ╚═════╝ ╚═╝  ╚═══╝
  //
  if (rerun_active)
  {

#define RERUN_DISPLAY 1
#if RERUN_DISPLAY
  std::string tcp_port = "127.0.0.1:9876";
// Initialize the rerun logger
#if RERUN_LOGGING < RERUN_OFF
  {
    rerun::RecordingStream recording{"FBGA", "FBGA"};
    recording.set_global();
  }
  rerun::RecordingStream &recording = rerun::RecordingStream::current();
  rerun::Error err = recording.connect_tcp(tcp_port, 10.0);
  err.throw_on_failure();
#endif
  // sleep 1s
  std::this_thread::sleep_for(std::chrono::seconds(1));
  //
  GG::integer numpoints = ceil(total_L);
  //
  std::vector<rerun::datatypes::Vec3D> points;
  for (GG::integer i = 0; i < numpoints; ++i) {
    GG::real ss = (GG::real) i * (total_L / (GG::real) (numpoints - 1));
    // Log the position of the trajectory
    points.emplace_back((GG::floating) traj_spline.eval("x_rl_m", ss),
                        (GG::floating) traj_spline.eval("y_rl_m", ss),
                        (GG::floating) traj_spline.eval("z_rl_m", ss));
  }
  //
  rerun::components::LineStrip3D line_strip(points);
  recording.log_static(
      "FBGA/track",
      rerun::LineStrips3D(line_strip)
          .with_colors(rerun::Color(colours["green"].r, colours["green"].g,
                                    colours["green"].b))
          .with_radii(rerun::Radius::scene_units(0.5)));
  //
  std::vector<rerun::Position3D> points_enve_up;
  std::vector<rerun::Position3D> points_enve_dw;
  points_enve_up.reserve(static_cast<std::size_t>(NUMPTSPLOTSHELL) * static_cast<std::size_t>(NUMPTSPLOTSHELL));
  points_enve_dw.reserve(static_cast<std::size_t>(NUMPTSPLOTSHELL) * static_cast<std::size_t>(NUMPTSPLOTSHELL));
  GG::NodeStruct3D fictitious_node;
  fictitious_node.s           = 0.0;
  fictitious_node.mu          = 0.0;
  fictitious_node.phi         = 0.0;
  fictitious_node.theta       = 0.0;
  fictitious_node.mu_prime    = 0.0;
  fictitious_node.phi_prime   = 0.0;
  fictitious_node.theta_prime = 0.0;
  fictitious_node.n           = 0.0;
  fictitious_node.chi         = 0.0;
  fbga_indy.eval_Omega_xyz_plot(fictitious_node);
  fbga_indy.eval_g_xyz_plot(fictitious_node);

  GG::output_plot_ggv_shell static_plot_shell = fbga_indy.eval_shell_plot( fictitious_node, 10.0,80.0,0.0);

  for(GG::size_type index = 0; index < static_plot_shell.v.size(); ++index) {
    points_enve_dw.emplace_back((GG::floating) static_plot_shell.a_tilde_y[index],
         (GG::floating) static_plot_shell.v[index],
         (GG::floating) static_plot_shell.a_tilde_x_min[index]);
    points_enve_up.emplace_back((GG::floating) static_plot_shell.a_tilde_y[index],
         (GG::floating) static_plot_shell.v[index],
         (GG::floating) static_plot_shell.a_tilde_x_max[index]);
  }

  recording.log_static(
      "FBGA/envelope/upper",
      rerun::Points3D(points_enve_up)
          .with_colors(rerun::Color(colours["red"].r, colours["red"].g,
                                    colours["red"].b)));
  recording.log_static(
      "FBGA/envelope/lower",
      rerun::Points3D(points_enve_dw)
          .with_colors(rerun::Color(colours["blue"].r, colours["blue"].g,
                                    colours["blue"].b)));
  //
  std::vector<rerun::Position3D> points_ggv_path;
  //
  points_ggv_path.reserve(static_cast<std::size_t>(numpts_gg));
  // 
  const GG::integer numpts_gg_plot = ceil(total_L*3.0);
  GG::real ds_plot = total_L / (GG::real) (numpts_gg_plot - 1);
  for (GG::integer i = 0; i < numpts_gg_plot; ++i) {
    GG::real space = (GG::real) i * ds_plot;
    //
    recording.set_time_seconds("space", space);
    //
    auto x_val = traj_spline.eval("x_rl_m", space);
    auto y_val = traj_spline.eval("y_rl_m", space);
    auto z_val = traj_spline.eval("z_rl_m", space);
    //
    auto vel_max_val = fbga_indy.eval_Vmax(space);
    auto Omega_x_val = fbga_indy.eval_Omega_x(space);
    auto Omega_y_val = fbga_indy.eval_Omega_y(space);
    auto Omega_z_val = fbga_indy.eval_Omega_z(space);
    //
    auto theta_spline_val = traj_spline.eval("theta_ref_rl_rad", space);
    auto phi_spline_val = traj_spline.eval("phi_ref_rl_rad", space);
    auto mu_spline_val = traj_spline.eval("mu_ref_rl_rad", space);
    auto dtheta_spline_val = traj_spline.eval("dtheta_ref_rl_radpm", space);
    auto dphi_spline_val = traj_spline.eval("dphi_ref_rl_radpm", space);
    auto dmu_spline_val = traj_spline.eval("dmu_ref_rl_radpm", space);
    auto chi_spline_val = traj_spline.eval("chi_rl_rad", space);
    auto n_spline_val = traj_spline.eval("n_rl_m", space);
    //
    auto segment_type_val = fbga_indy.eval_segment_type(space);
    //
    auto v_val = fbga_indy.eval_V(space);
    auto a_hat_x_val = fbga_indy.eval_A_hat_x(space);
    auto a_hat_y_val = fbga_indy.eval_A_hat_y(space);
    auto a_hat_z_val = fbga_indy.eval_A_hat_z(space);
    //
    auto a_tilde_x_val = fbga_indy.eval_A_tilde_x(space);
    auto a_tilde_y_val = fbga_indy.eval_A_tilde_y(space);
    auto a_tilde_z_val = fbga_indy.eval_A_tilde_z(space);
    //
    auto v_dot_val = fbga_indy.eval_V_dot(space);
    //
    auto v_mpc_val = traj_spline.eval("v_rl_mps", space);
    //
    auto a_x_mpc_val = traj_spline.eval("ax_rl_mps2", space);
    auto a_y_mpc_val = traj_spline.eval("ay_rl_mps2", space);
    //
    auto alpha_val = fbga_indy.eval_alpha(space);
    // Log values
    recording.log("x", rerun::Scalar{x_val});
    recording.log("y", rerun::Scalar{y_val});
    recording.log("z", rerun::Scalar{z_val});
    //
    recording.log("vel_max", rerun::Scalar{vel_max_val});
    recording.log("Omega_x", rerun::Scalar{Omega_x_val});
    recording.log("Omega_y", rerun::Scalar{Omega_y_val});
    recording.log("Omega_z", rerun::Scalar{Omega_z_val});
    //
    recording.log("theta_spline", rerun::Scalar{theta_spline_val});
    recording.log("phi_spline", rerun::Scalar{phi_spline_val});
    recording.log("mu_spline", rerun::Scalar{mu_spline_val});
    recording.log("dtheta_spline", rerun::Scalar{dtheta_spline_val});
    recording.log("dphi_spline", rerun::Scalar{dphi_spline_val});
    recording.log("dmu_spline", rerun::Scalar{dmu_spline_val});
    recording.log("chi_spline", rerun::Scalar{chi_spline_val});
    recording.log("n_spline", rerun::Scalar{n_spline_val});
    //
    recording.log("v", rerun::Scalar{v_val});
    recording.log("a_hat_x", rerun::Scalar{a_hat_x_val});
    recording.log("a_hat_y", rerun::Scalar{a_hat_y_val});
    recording.log("a_hat_z", rerun::Scalar{a_hat_z_val});
    //
    recording.log("a_tilde_x", rerun::Scalar{a_tilde_x_val});
    recording.log("a_tilde_y", rerun::Scalar{a_tilde_y_val});
    recording.log("a_tilde_z", rerun::Scalar{a_tilde_z_val});
    //
    recording.log("v_dot", rerun::Scalar{v_dot_val});
    //
    recording.log("segment_type", rerun::Scalar{static_cast<int>(segment_type_val)});
    //
    recording.log("v_mpc", rerun::Scalar{v_mpc_val});
    recording.log("a_x_mpc", rerun::Scalar{a_x_mpc_val});
    recording.log("a_y_mpc", rerun::Scalar{a_y_mpc_val});
    //
    recording.log("alpha", rerun::Scalar{alpha_val});
    //
    // Only execute every 5 iterations to reduce computation and logging
    if (i % 1 == 0) {
      fictitious_node.s           = space;
      fictitious_node.mu          = mu_spline_val;
      fictitious_node.phi         = phi_spline_val;
      fictitious_node.theta       = theta_spline_val;
      fictitious_node.mu_prime    = dmu_spline_val;
      fictitious_node.phi_prime   = dphi_spline_val;
      fictitious_node.theta_prime = dtheta_spline_val;
      fictitious_node.n           = n_spline_val;
      fictitious_node.chi         = chi_spline_val;
      fictitious_node.alpha       = alpha_val;
      fbga_indy.eval_Omega_xyz_plot(fictitious_node);
      fbga_indy.eval_g_xyz_plot(fictitious_node);
      //
      GG::output_plot_ggv_shell dynamic_plot_shell = fbga_indy.eval_shell_plot( fictitious_node, 10.0,80.0,0.0);
      //
      points_enve_dw.clear();
      points_enve_up.clear();
      points_enve_dw.reserve(static_cast<std::size_t>(NUMPTSPLOTSHELL) * static_cast<std::size_t>(NUMPTSPLOTSHELL));
      points_enve_up.reserve(static_cast<std::size_t>(NUMPTSPLOTSHELL) * static_cast<std::size_t>(NUMPTSPLOTSHELL));
      for(GG::size_type index = 0; index < dynamic_plot_shell.v.size(); ++index) {
      points_enve_dw.emplace_back((GG::floating) dynamic_plot_shell.a_tilde_y[index],
        (GG::floating) dynamic_plot_shell.v[index],
        (GG::floating) dynamic_plot_shell.a_tilde_x_min[index]);
      points_enve_up.emplace_back((GG::floating) dynamic_plot_shell.a_tilde_y[index],
        (GG::floating) dynamic_plot_shell.v[index],
        (GG::floating) dynamic_plot_shell.a_tilde_x_max[index]);
      }
      //
      recording.log(
        "FBGA/envelope/upper_moving",
        rerun::Points3D(points_enve_up)
        .with_colors(rerun::Color(colours["mint"].r, colours["mint"].g, colours["mint"].b)));
      recording.log(
        "FBGA/envelope/lower_moving",
        rerun::Points3D(points_enve_dw)
          .with_colors(rerun::Color(colours["mint"].r, colours["mint"].g, colours["mint"].b)));
    }
    //
    recording.log(
      "FBGA/GGv_point",
      rerun::Points3D( rerun::Position3D( (GG::floating)a_tilde_y_val, (GG::floating)v_val, (GG::floating)a_tilde_x_val) )
          .with_colors(rerun::Color(colours["cyan"].r, colours["cyan"].g,colours["cyan"].b))
          .with_radii(0.5f));
    //
    recording.log(
      "FBGA/Traj_point",
      rerun::Points3D( rerun::Position3D( (GG::floating)x_val, (GG::floating)y_val, (GG::floating)z_val) )
          .with_colors(rerun::Color(colours["cyan"].r, colours["cyan"].g,colours["cyan"].b))
          .with_radii(3.5f));
    //
    points_ggv_path.emplace_back((GG::floating)a_tilde_y_val, (GG::floating)v_val, (GG::floating)a_tilde_x_val);
  }
  //
  recording.log_static(
      "FBGA/GGv_path",
      rerun::Points3D(points_ggv_path)
          .with_colors(rerun::Color(colours["green"].r, colours["green"].g,
                                    colours["green"].b)));
  //
  std::this_thread::sleep_for(std::chrono::seconds(1));
#endif
  //

  } // end rerun_active


  return 0;
}
