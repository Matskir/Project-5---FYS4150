#include <iostream>
#include <complex>
#include <armadillo>


//Needed to use 1i as the complex i
using namespace std::complex_literals;

//Function that calculates the index k from the indeces i and j
int ind_k(int ind_i, int ind_j, int M){
  return ind_i + (ind_j-1)*(M-2) -1;
}


//Function that fills either the matrix A or B with values according to Crank-Nicolson
void mat_A_B(double M, arma::cx_double r, arma::cx_vec a, arma::sp_cx_mat& matrix){
  matrix.diag(0) = a;
  matrix.diag(1).fill(-r);
  matrix.diag(-1).fill(-r);
  matrix.diag(M-2).fill(-r);
  matrix.diag(-(M-2)).fill(-r);

  for (int i=0; i<= matrix.diag(1).n_rows; i++){

    if (i%(int)(M-2) == 0 && i!= 0){
      matrix.diag(1)(i-1) = 0.;
      matrix.diag(-1)(i-1) = 0.;
    }
  }
}


//Function that fill both A and B according to Crank-Nicolson
void fill_matrices(double M, double h, double dt, arma::sp_mat V, arma::sp_cx_mat& A, arma::sp_cx_mat& B){
  arma::cx_vec a = arma::cx_vec((M-2)*(M-2));
  arma::cx_vec b = arma::cx_vec((M-2)*(M-2));
  arma::cx_double r = 1i*dt/(2.*h*h);


  for (int i=1; i < (M-2); i++){
    for (int j=1; j < (M-2); j++){

      int k = ind_k(i,j,M);

      a(k) = 1. + 4.*r + 1i*dt*((double) V(i,j))/2.;
      b(k) = 1. - 4.*r - 1i*dt*((double) V(i,j))/2.;
    }
  }
  mat_A_B(M,r,a,A);
  mat_A_B(M,-r,b,B);
}


//Function that computes the "wave function" at n+1
arma::cx_vec next_u(arma::sp_cx_mat A, arma::sp_cx_mat B, arma::cx_vec u){
  arma::cx_vec b = B*u;
  arma::cx_vec u_next = arma::spsolve(A,b);
  return u_next;
}

//Function that sets up the initial state u0
arma::cx_mat set_up_initial_state(double xc, double yc, double sigma_x, double sigma_y, double px, double py, double M, double h){
  arma::cx_mat u_in = arma::cx_mat(M,M);

  for (int i=0; i<=M-1; i++){
    for (int j=0; j<=M-1; j++){
      u_in(j,i) = exp(-(i*h-xc)*(i*h-xc)/(2*sigma_x*sigma_x) - (j*h-yc)*(j*h-yc)/(2*sigma_y*sigma_y) + 1i*px*(i*h-xc) + 1i*py*(j*h-yc));
    }
  }

  //boundary conditions
  int n = u_in.n_cols;
  u_in.col(0) = u_in.col(0)*0.;
  u_in.col(n-1) = u_in.col(n-1)*0.;

  int n_row = u_in.n_rows;
  u_in.row(0) = u_in.row(0)*0;
  u_in.row(n-1) = u_in.row(n-1)*0;
return u_in/std::sqrt(arma::sum(u_in.as_col()%arma::conj(u_in.as_col())).real());
}


//Function that sets up the double slit potential
arma::sp_mat set_up_potential(double v0, double h, double M){
  int ind_thickness_x = 0.02/h;
  int ind_wall_pos_x_left = 0.5/h-0.5*ind_thickness_x;
  int ind_y_centre = 0.5/h;
  int ind_separation_y = 0.05/h;
  int ind_slit_ap_y = 0.05/h;

  arma::sp_mat V = arma::sp_mat((M-2),(M-2));

  for (int i=0; i<M-2; i++){
    for (int j=0; j<M-2; j++){
      if (ind_wall_pos_x_left <= i && i <= ind_wall_pos_x_left+ind_thickness_x){
        V(j,i) = v0;

        if (ind_y_centre-0.5*ind_separation_y-ind_slit_ap_y <= j && j <= ind_y_centre-0.5*ind_separation_y){
          V(j,i) = 0;
        }
        else if (ind_y_centre+0.5*ind_separation_y <= j && j<= ind_y_centre+0.5*ind_separation_y+ind_slit_ap_y){
          V(j,i) = 0;
        }
      }
    }
  }
  return V;
}

//Function that sets up the single slit potential
arma::sp_mat set_up_single_slit(double v0, double h, double M){
  int ind_thickness_x = 0.02/h;
  int ind_wall_pos_x_left = 0.5/h-0.5*ind_thickness_x;
  int ind_y_centre = 0.5/h;
  int ind_separation_y = 0.05/h;
  int ind_slit_ap_y = 0.05/h;

  arma::sp_mat V = arma::sp_mat((M-2),(M-2));

  for (int i=0; i<M-2; i++){
    for (int j=0; j<M-2; j++){
      if (ind_wall_pos_x_left <= i && i <= ind_wall_pos_x_left+ind_thickness_x){
        V(j,i) = v0;

        if (ind_y_centre-0.5*ind_slit_ap_y <= j && j <= ind_y_centre+0.5*ind_slit_ap_y){
          V(j,i) = 0;
        }
      }
    }
  }
  return V;
}

//Function that sets up the triple slit potential
arma::sp_mat set_up_triple_slit(double v0, double h, double M){
  int ind_thickness_x = 0.02/h;
  int ind_wall_pos_x_left = 0.5/h-0.5*ind_thickness_x;
  int ind_y_centre = 0.5/h;
  int ind_separation_y = 0.05/h;
  int ind_slit_ap_y = 0.05/h;

  arma::sp_mat V = arma::sp_mat((M-2),(M-2));

  for (int i=0; i<M-2; i++){
    for (int j=0; j<M-2; j++){
      if (ind_wall_pos_x_left <= i && i <= ind_wall_pos_x_left+ind_thickness_x){
        V(j,i) = v0;

        if (ind_y_centre-0.5*ind_slit_ap_y <= j && j <= ind_y_centre+0.5*ind_slit_ap_y){
          V(j,i) = 0;
        }

        if (ind_y_centre+0.5*ind_slit_ap_y+ind_separation_y <= j && j <= ind_y_centre+1.5*ind_slit_ap_y+ind_separation_y){
          V(j,i) = 0;
        }

        if (ind_y_centre-1.5*ind_slit_ap_y-ind_separation_y <= j && j <= ind_y_centre-0.5*ind_slit_ap_y-ind_separation_y){
          V(j,i) = 0;
        }
      }
    }
  }
  return V;

}
