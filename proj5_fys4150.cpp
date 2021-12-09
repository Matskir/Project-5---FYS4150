#include "header_proj5.hpp"

int main(){
  //Declare and define all parameters
  double h = 0.005;
  double dt = 2.5e-5;
  double T = 0.002;
  double xc = 0.25;
  double sigma_x = 0.05;
  double px = 200;
  double yc = 0.5;
  double sigma_y = 0.2;
  double py = 0;
  double v0 = 1e10;
  double M = 1./h;


  //No potential
  //arma::sp_mat V = arma::sp_mat((M-2),(M-2));

  //Double slit potential
  arma::sp_mat V = set_up_potential(v0, h, M);

  //Single slit potential
  //arma::sp_mat V = set_up_single_slit(v0, h, M);

  //Triple slit potential
  //arma::sp_mat V = set_up_triple_slit(v0, h, M);

  //Declare matrices A and B
  arma::sp_cx_mat A = arma::sp_cx_mat((M-2)*(M-2),(M-2)*(M-2));
  arma::sp_cx_mat B = arma::sp_cx_mat((M-2)*(M-2),(M-2)*(M-2));

  //Fill matrices A and B
  fill_matrices(M,h,dt,V,A,B);

  //Compute the initial state u0 according to defined parameters
  arma::cx_mat u_init = set_up_initial_state(xc, yc, sigma_x, sigma_y, px, py, M, h);

  //Remove boundaries so as to only include internal points
  arma::uvec index_row = {0,M-1};
  arma::uvec index_col = {0,M-1};
  u_init.shed_cols(index_row);
  u_init.shed_rows(index_col);

  //Declare a data cube to be filled with the wave functions at each time n
  arma::cx_cube time_cube = arma::cx_cube(M-2,M-2,T/dt);
  time_cube.slice(0) = u_init;

  //Compute u at n+1 and add it to data cube
  for (int i=1; i<T/dt; i++){
    arma::cx_mat u_prev = time_cube.slice(i-1);
    arma::cx_vec u_next =  next_u(A,B,u_prev.as_col());
    arma::cx_mat u_next_mat = reshape(u_next,M-2,M-2);
    time_cube.slice(i) = u_next_mat;
  }
  //save data cube
  time_cube.save("cube.bin");
  return 0;
}
