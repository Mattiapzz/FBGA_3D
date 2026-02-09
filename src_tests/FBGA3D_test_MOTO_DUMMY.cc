#include <iostream>

#include "FBGA3D/gggv_MOTO.hxx"
#include "FBGA3D/gg_utils.hxx"
#include "FBGA3D/types.hxx"


// chrono


int main()

{

  // 

  std::cout << "FBGA3D MOTO test started\n";
  std::cout << "Example class GGGV MOTO\n";


  //  ██████╗  ██████╗  ██████╗ ██╗   ██╗
  // ██╔════╝ ██╔════╝ ██╔════╝ ██║   ██║
  // ██║  ███╗██║  ███╗██║  ███╗██║   ██║
  // ██║   ██║██║   ██║██║   ██║╚██╗ ██╔╝
  // ╚██████╔╝╚██████╔╝╚██████╔╝ ╚████╔╝ 
  //  ╚═════╝  ╚═════╝  ╚═════╝   ╚═══╝  
                                      

  GG::gggv_MOTO gggv_moto;

  //
  std::cout << "======================== FBGA STATS ========================\n";
  std::cout << "=============================================================\n";

  // ███╗   ██╗██╗   ██╗███╗   ███╗███████╗██████╗ ██╗ ██████╗
  // ████╗  ██║██║   ██║████╗ ████║██╔════╝██╔══██╗██║██╔════╝
  // ██╔██╗ ██║██║   ██║██╔████╔██║█████╗  ██████╔╝██║██║     
  // ██║╚██╗██║██║   ██║██║╚██╔╝██║██╔══╝  ██╔══██╗██║██║     
  // ██║ ╚████║╚██████╔╝██║ ╚═╝ ██║███████╗██║  ██║██║╚██████╗
  // ╚═╝  ╚═══╝ ╚═════╝ ╚═╝     ╚═╝╚══════╝╚═╝  ╚═╝╚═╝ ╚═════╝
                                                            

  std::cout << "gggv_MOTO object created\n";

  // evaluate engine for different velocities
  std::cout << "Evaluating engine for different velocities:\n";
  for (GG::real V = 10.0; V <= 100.0; V += 5.0) {
    GG::real a_x_eng = gggv_moto.a_x_eng(V);
    std::cout << "V: " << V << " m/s, a_x_eng: " << a_x_eng << " m/s^2\n";
  }

  std::cout << "=============================================================\n";

  // fix velocity 10 m/s eval a_x_max and a_x_min for values between -a_y_lim and a_y_lim
  GG::real V = 10.0;
  GG::real az_tilde = 9.81; // assuming az_tilde is g
  std::cout << "Evaluating a_x_max and a_x_min for V = " << V << " m/s and az_tilde = " << az_tilde << " m/s^2:\n";
  GG::real a_y_lim = gggv_moto.a_y_lim(V, az_tilde);
  for (GG::real a_y_tilde = -a_y_lim; a_y_tilde <= a_y_lim; a_y_tilde += 0.5) {
    GG::real a_x_max = gggv_moto.a_x_push(a_y_tilde, V, az_tilde);
    GG::real a_x_min = gggv_moto.a_x_pull(a_y_tilde, V, az_tilde);
    std::cout << "a_y_tilde: " << a_y_tilde << " m/s^2, a_x_max: " << a_x_max << " m/s^2, a_x_min: " << a_x_min << " m/s^2\n";
  }
  std::cout << "=============================================================\n";
  //
  az_tilde = 8.0; // changing az_tilde to 8.0 m/s^2
  std::cout << "Evaluating a_x_max and a_x_min for V = " << V << " m/s and az_tilde = " << az_tilde << " m/s^2:\n";
  a_y_lim = gggv_moto.a_y_lim(V, az_tilde);
  for (GG::real a_y_tilde = -a_y_lim; a_y_tilde <= a_y_lim; a_y_tilde += 0.5) {
    GG::real a_x_max = gggv_moto.a_x_push(a_y_tilde, V, az_tilde);
    GG::real a_x_min = gggv_moto.a_x_pull(a_y_tilde, V, az_tilde);
    std::cout << "a_y_tilde: " << a_y_tilde << " m/s^2, a_x_max: " << a_x_max << " m/s^2, a_x_min: " << a_x_min << " m/s^2\n";
  }
  std::cout << "=============================================================\n";

  //
  az_tilde = 11.0; 
  std::cout << "Evaluating a_x_max and a_x_min for V = " << V << " m/s and az_tilde = " << az_tilde << " m/s^2:\n";
  a_y_lim = gggv_moto.a_y_lim(V, az_tilde);
  for (GG::real a_y_tilde = -a_y_lim; a_y_tilde <= a_y_lim; a_y_tilde += 0.5) {
    GG::real a_x_max = gggv_moto.a_x_push(a_y_tilde, V, az_tilde);
    GG::real a_x_min = gggv_moto.a_x_pull(a_y_tilde, V, az_tilde);
    std::cout << "a_y_tilde: " << a_y_tilde << " m/s^2, a_x_max: " << a_x_max << " m/s^2, a_x_min: " << a_x_min << " m/s^2\n";
  }
  std::cout << "=============================================================\n";






  return 0;
}
